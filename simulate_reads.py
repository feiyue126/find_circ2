#!/usr/bin/env python

__version__ = "0.90"
__author__ = "Marvin Jens"
__credits__ = ["Marvin Jens","Marcel Schilling","Petar Glazar","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marvin.jens@mdc-berlin.de"

from byo.gene_model import transcripts_from_UCSC, CircRNA
from byo import rev_comp
import os,sys
import optparse
import numpy as np
import logging
from collections import defaultdict
from gzip import GzipFile
lin_counts = defaultdict(int)
circ_counts = defaultdict(int)

usage = """
   cat circ_models.ucsc | %prog [options] > interlaced_paired_end_reads.fa
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome FASTA file")
parser.add_option("","--mutate",dest="mut_rate",type=float,default=0,help="per base mutation rate, between 0 and 1 (default=0)")
parser.add_option("","--n-frags",dest="n_frags",type=int,default=100,help="number of fragments to simulate (default=100)")
parser.add_option("","--fpk",dest="fpk",action="store_true",default=False, help="if set, --n-frags is interpreted as frags-per-kilobase")
parser.add_option("","--frag-len",dest="frag_len",type=int,default=350,help="fragment length to simulate (default=350)")
parser.add_option("","--read-len",dest="read_len",type=int,default=100,help="read length to simulate (default=100)")
parser.add_option("-o","--output",dest="output",default="simulation",help="path, where to store the output (default='./simulation')")
parser.add_option("","--stdout",dest="stdout",default=None,choices=['circs','lins','reads','multi','test'],help="use to direct chosen type of output (circs, lins, reads, multi) to stdout instead of file")

options,args = parser.parse_args()

if not (options.genome or options.system):
    error("need to specify either model system database (-S) or genome FASTA file (-G).")
    sys.exit(1)

# prepare output files
if not os.path.isdir(options.output):
    os.mkdir(options.output)

# prepare logging system
FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
logging.basicConfig(level=logging.INFO,format=FORMAT,filename=os.path.join(options.output,"run.log"),filemode='w')
logger = logging.getLogger('simulate_reads.py')
logger.info("simulate_reads.py {0} invoked as '{1}'".format(__version__," ".join(sys.argv)))

circs_file = file(os.path.join(options.output,"circ_splice_sites.bed"),"w")
lins_file  = file(os.path.join(options.output,"lin_splice_sites.bed"),"w")
reads_file = GzipFile(os.path.join(options.output,"simulated_reads.fa.gz"),"w")

# redirect output to stdout, if requested
if options.stdout:
    varname = "{0}_file".format(options.stdout)
    out_file = globals()[varname]
    out_file.write('# redirected to stdout\n')
    logger.info('redirected {0} to stdout'.format(options.stdout))
    globals()[varname] = sys.stdout

if options.system:
    import importlib
    system = importlib.import_module("byo.systems.{options.system}".format(**locals()))
else:
    from byo.track import load_track
    class S(object):
        pass
    
    system = S()
    system.genome = load_track(options.genome)

def store_splices(data, dst, prefix="sim"):
    for i,(coord,count) in enumerate(sorted(data.items())):
        chrom, start, end, strand = coord
        name = "{0}_{1}".format(prefix,i)
        out = [chrom, start, end, name, count, strand]
        
        dst.write('\t'.join([str(o) for o in out]) + '\n')
    
def test_str(mate, rec_lin = {}, rec_circ = {}):
    ori = mate.origin
    seg_bounds = np.roll(mate.exon_bounds, -mate.map_to_exon(mate.origin_spliced+1), axis=0)

    parts = []
    last_s, last_e = seg_bounds[0]
    for i,(s,e) in enumerate(seg_bounds[::mate.dir]):
        if i == 0:
            # first segment: yield origin
            parts.append("O:{chrom}:{ori}:{sense};M:{m}".format(chrom=mate.chrom, ori=mate.origin, sense=mate.sense, m=e-s) )
        elif s > last_s:
            # we move forward: linear splicing
            parts.append("LS:{start}:{end};M:{m}".format(start=last_e - ori, end=s - ori, m=e-s) )
            rec_lin[ (mate.chrom, last_e, s, mate.sense) ] += 1

        elif s < last_s:
            # we move backward: circRNA backsplicing
            parts.append("CS:{start}:{end};M:{m}".format(start=s - ori, end=last_e - ori, m=e-s) )
            rec_circ[ (mate.chrom, s, last_e, mate.sense) ] += 1
    
        last_s, last_e = s,e

    return ";".join(parts)

def mutate(seq, rate):
    other = {
        'a' : ['C','G','T'],
        'c' : ['A','G','T'],
        'g' : ['A','C','T'],
        't' : ['A','C','G'],
    }
    
    seq = list(seq.lower())
    p = np.random.rand(len(seq))
    for i in (p < rate).nonzero()[0]:
        seq[i] = other[seq[i]][np.random.randint(0,3)]
    
    return "".join(seq)

for circ in transcripts_from_UCSC(sys.stdin, system=system, tx_type=CircRNA):

    L = circ.spliced_length
    
    if options.fpk:
        n = int(options.n_frags * L / 1E3)
    else:
        n = options.n_frags

    logger.info("simulating {0} fragments for {1}".format(n, circ.name))

    rnd_pos = np.random.randint(0, L, n)
    for i,frag_start in enumerate(rnd_pos):
        
        m1_start, m1_end = frag_start, frag_start + options.read_len
        m2_start, m2_end = frag_start + options.frag_len - options.read_len, frag_start + options.frag_len
        
        m1_g_start, m1_g_end = circ.map_from_spliced(m1_start), circ.map_from_spliced(m1_end)
        m2_g_start, m2_g_end = circ.map_from_spliced(m2_start), circ.map_from_spliced(m2_end)
        
        mate1 = circ.cut(m1_g_start, m1_g_end)
        mate2 = circ.cut(m2_g_start, m2_g_end)
        
        assert mate1.spliced_length == options.read_len
        assert mate2.spliced_length == options.read_len
        
        # simulate forward/reverse mate pairs
        mate1_seq = mate1.spliced_sequence
        mate2_seq = rev_comp(mate2.spliced_sequence)
        
        if options.mut_rate:
            mate1_seq = mutate(mate1_seq, options.mut_rate)
            mate2_seq = mutate(mate2_seq, options.mut_rate)
            
        m1_str = test_str(mate1, rec_lin = lin_counts, rec_circ = circ_counts)
        m2_str = test_str(mate2, rec_lin = lin_counts, rec_circ = circ_counts)

        read_name = "sim_{circ.name}_{i}___{m1_str}|{m2_str}".format(**locals())
        reads_file.write(">{read_name}\n{mate1_seq}\n>{read_name}\n{mate2_seq}\n".format(**locals()) )

logger.info("storing {0} linear junction span counts".format(len(lin_counts)))
store_splices(lin_counts, lins_file, "sim_lin")

logger.info("storing {0} circular junction span counts".format(len(circ_counts)))
store_splices(circ_counts, circs_file, "sim_circ")
