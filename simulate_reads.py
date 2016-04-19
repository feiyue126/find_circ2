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

# store the number of [n_weighted, n_frag, n_span, min_span]
lin_counts = defaultdict( lambda : np.zeros(4) )
circ_counts = defaultdict( lambda : np.zeros(4) )

usage = """
   cat circ_models.ucsc | %prog [options] > interlaced_paired_end_reads.fa
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome FASTA file")
parser.add_option("","--mutate",dest="mut_rate",type=float,default=0,help="per base mutation rate, between 0 and 1 (default=0)")
parser.add_option("","--min-segment",dest="min_seg",type=int,default=13,help="minimum size of span to be detectable/~BWA MEM min. segment size (default=13)")
parser.add_option("","--n-frags",dest="n_frags",type=int,default=100,help="number of fragments to simulate (default=100)")
parser.add_option("","--fpk",dest="fpk",action="store_true",default=False, help="if set, --n-frags is interpreted as frags-per-kilobase")
parser.add_option("","--frag-len",dest="frag_len",type=int,default=350,help="fragment length to simulate (default=350)")
parser.add_option("","--read-len",dest="read_len",type=int,default=100,help="read length to simulate (default=100)")
parser.add_option("","--seed",dest="seed",type=int,default=0,help="seed for pseudo random number generator (default=random)")
parser.add_option("-o","--output",dest="output",default="simulation",help="path, where to store the output (default='./simulation')")
parser.add_option("","--stdout",dest="stdout",default=None,choices=['circs','lins','reads','multi','test'],help="use to direct chosen type of output (circs, lins, reads, multi) to stdout instead of file")

options,args = parser.parse_args()

# prepare output directory
if not os.path.isdir(options.output):
    os.makedirs(options.output)

# prepare logging system
FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
logging.basicConfig(level=logging.INFO,format=FORMAT,filename=os.path.join(options.output,"run.log"),filemode='w')
logger = logging.getLogger('simulate_reads.py')
logger.info("simulate_reads.py {0} invoked as '{1}'".format(__version__," ".join(sys.argv)))

if not (options.genome or options.system):
    print "need to specify either model system database (-S) or genome FASTA file (-G)."
    sys.exit(1)

# prepare output files
circs_file = file(os.path.join(options.output,"circ_splice_sites.bed"),"w")
lins_file  = file(os.path.join(options.output,"lin_splice_sites.bed"),"w")
reads_file = GzipFile(os.path.join(options.output,"simulated_reads.fa.gz"),"w")

## ugly global variables to hold differnt parts together
# the name dictionary allows the BED-file writer in the 
# end to use the correct name for each junction
circ_junction_names = {}

# redirect output to stdout, if requested
if options.stdout:
    varname = "{0}_file".format(options.stdout)
    out_file = globals()[varname]
    out_file.write('# redirected to stdout\n')
    logger.info('redirected {0} to stdout'.format(options.stdout))
    globals()[varname] = sys.stdout

# initialize random number generator
if options.seed:
    np.random.seed(options.seed)

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
    dst.write('\t'.join(["#chrom","start","end","name","n_detectable","strand","n_span","n_weight","n_frags"]) + '\n')
    for i,(coord,(n_weight, n_frags, n_span, n_detect)) in enumerate(sorted(data.items())):
        chrom, start, end, strand = coord
        # get circ-name if this a circ junction, or generate generic name
        name = circ_junction_names.get(coord, "{0}_{1}".format(prefix,i))
        out = [chrom, start, end, name, n_detect, strand, n_span, n_weight, n_frags]
        
        dst.write('\t'.join([str(o) for o in out]) + '\n')
    
def test_str(mate, lin_juncs = {}, circ_juncs = {}, circ=None):
    splices = 0
    lin_detectable = defaultdict(int)
    parts = []
    lin_span = defaultdict(int)

    last_s, last_e = mate.start, mate.end
    ori = mate.start
    # always represent in + orientation of the chromosome!
    for i,(s,e) in enumerate(mate.exon_bounds):
        if i == 0:
            # first segment: yield origin
            parts.append("O:{chrom}:{ori}:{sense};M:{m}".format(chrom=mate.chrom, ori=ori, sense=mate.sense, m=e-s) )
        else:
            # we must have gotten here by splicing
            if s == circ.start or e == circ.end:
                # these are circRNA bounds, not linked by linear splicing: skip for now!
                pass
            else:
                parts.append("LS:{start}:{end};M:{m}".format(start=last_e - ori, end=s - ori, m=e-s) )
                coord = (mate.chrom, last_e, s, mate.sense)

                left, right = mate.exon_lengths[i-1:i+1]
                # record the number of nucleotides that span the junction
                span = min(left, right)
                lin_span[coord] = max(lin_span[coord], span)
   
        last_s, last_e = s,e

    for coord, span in lin_span.items():
        if span >= options.min_seg:
            lin_juncs[coord] += np.array([1,1])
        else:
            lin_juncs[coord] += np.array([1,0])
            
    if mate.start == circ.start and mate.end == circ.end:
        parts.append("CS:{start}:{end};M:{m}".format(start=circ.start - ori, end=circ.end - ori, m=e-s) )
        span = min(mate.exon_lengths[0], mate.exon_lengths[-1])
        coord = (mate.chrom, circ.start, circ.end, mate.sense)
        if span >= options.min_seg:
            circ_juncs[coord] += np.array([1,1]) 
        else:
            lin_juncs[coord] += np.array([1,0])
            
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

def test_mate(mate, circ):
    """
    tests if the exon boundaries of the mate make sense
    """
    circseq = circ.spliced_sequence.lower() * 2
    res = mate.spliced_sequence.lower() in circseq
    if not res:
        print "WTF?"
        print "mate",mate
        print mate.c_start_end
        s,e = mate.c_start_end
        print e-s, circ.spliced_length
        print "circ",circ
        print circseq
        print mate.spliced_sequence.lower()
        print res

    return res

for circ in transcripts_from_UCSC(sys.stdin, system=system, tx_type=CircRNA):

    # keep circ name by junction coordinate for later, when we write reads and need the names!
    circ_junction_names[ (circ.chrom, circ.start, circ.end, circ.sense) ] = circ.name

    L = circ.spliced_length
    if L <= options.read_len:
        logger.warning("skipping {0} because it is shorter than read length, which is currently not supported".format(circ.name))
        continue

    # TODO: only for now to make sure that circs are comparable between different red-length runs
    if L <= 200:
        logger.warning("skipping {0} because it is shorter than read length, which is currently not supported".format(circ.name))
        continue

    #print circ.name, L
    
    if options.fpk:
        n = int(options.n_frags * L / 1E3)
    else:
        n = options.n_frags

    logger.info("simulating {0} fragments for {1}".format(n, circ.name))

    rnd_pos = np.random.randint(0, L, n)
    #rnd_pos = [0,100,200,400]
    for i,frag_start in enumerate(rnd_pos):
        #frag_start = circ.map_to_spliced(16323340)
        #print ">>> frag_start",frag_start
        m1_start, m1_end = frag_start, frag_start + options.read_len
        m2_start, m2_end = frag_start + options.frag_len - options.read_len, frag_start + options.frag_len
        
        m1_g_start, m1_g_end = circ.map_from_spliced(m1_start), circ.map_from_spliced(m1_end)
        m2_g_start, m2_g_end = circ.map_from_spliced(m2_start), circ.map_from_spliced(m2_end)
        
        #print circ.name, L
        #print "M1", m1_start, m1_end, m1_g_start, m1_g_end
        mate1 = circ.cut(m1_g_start, m1_g_end)
        mate1.c_start_end = m1_start, m1_end
        #print "M2", m2_start, m2_end, m2_g_start, m2_g_end
        mate2 = circ.cut(m2_g_start, m2_g_end)
        mate2.c_start_end = m2_start, m2_end
        
        assert mate1.spliced_length == options.read_len
        assert mate2.spliced_length == options.read_len
        assert test_mate(mate1, circ)
        assert test_mate(mate2, circ)
        #if not test_mate(mate2, circ):
            #print mate2
            #print "origin",mate2.origin, mate2.origin_spliced
            #print mate2.spliced_sequence
            #print circ
            #sys.exit(1)
        
        # simulate forward/reverse mate pairs
        mate1_seq = mate1.spliced_sequence
        mate2_seq = rev_comp(mate2.spliced_sequence)
        
        if options.mut_rate:
            mate1_seq = mutate(mate1_seq, options.mut_rate)
            mate2_seq = mutate(mate2_seq, options.mut_rate)
            
        lin_juncs = defaultdict(lambda : np.zeros(2) )
        circ_juncs = defaultdict(lambda : np.zeros(2) )
            
        m1_str = test_str(mate1, lin_juncs, circ_juncs, circ=circ)
        m2_str = test_str(mate2, lin_juncs, circ_juncs, circ=circ)

        splices = len(lin_juncs) + len(circ_juncs)
        if splices > 1:
            w = 1./(splices-1)
        else:
            w = 1
            
        # store the counts
        for coord, (counts, detectable) in lin_juncs.items():
            lin_counts[coord] += np.array([w, 1, counts, detectable > 0])

        is_detectable_circ = False
        for coord, (counts, detectable) in circ_juncs.items():
            if detectable:
                is_detectable_circ = True
            circ_counts[coord] += np.array([w, 1, counts, detectable > 0])
        
        if is_detectable_circ:
            flags = "_DETECTABLE"
        else:
            flags = ""

        #if circ_detect_long:
            #parts.append('LONG_DETECT')

        read_name = "{circ.name}_sim_{i}{flags}___{m1_str}|{m2_str}".format(**locals())
        reads_file.write(">{read_name}\n{mate1_seq}\n>{read_name}\n{mate2_seq}\n".format(**locals()) )

logger.info("storing {0} linear junction span counts".format(len(lin_counts)))
store_splices(lin_counts, lins_file, "lin_sim")

logger.info("storing {0} circular junction span counts".format(len(circ_counts)))
store_splices(circ_counts, circs_file, "circ_sim")
