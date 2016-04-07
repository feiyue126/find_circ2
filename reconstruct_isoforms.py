#!/usr/bin/env python

__version__ = "0.90"
__author__ = "Marvin Jens"
__credits__ = ["Marvin Jens","Marcel Schilling","Petar Glazar","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marvin.jens@mdc-berlin.de"

from byo.gene_model import transcripts_from_UCSC, CircRNA, Transcript
from byo import rev_comp
import os,sys
import optparse
import numpy as np
import bisect
import logging
from collections import defaultdict
from gzip import GzipFile

usage = """
   %prog [options] multi_events.tsv
"""

parser = optparse.OptionParser(usage=usage)
#parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
#parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome FASTA file")
parser.add_option("","--known-exons",dest="known_exons",default="",help="GTF file with known exons")
#parser.add_option("","--min-segment",dest="min_seg",type=int,default=13,help="minimum size of span to be detectable/~BWA MEM min. segment size (default=13)")
#parser.add_option("","--n-frags",dest="n_frags",type=int,default=100,help="number of fragments to simulate (default=100)")
#parser.add_option("","--fpk",dest="fpk",action="store_true",default=False, help="if set, --n-frags is interpreted as frags-per-kilobase")
#parser.add_option("","--frag-len",dest="frag_len",type=int,default=350,help="fragment length to simulate (default=350)")
#parser.add_option("","--read-len",dest="read_len",type=int,default=100,help="read length to simulate (default=100)")
parser.add_option("-o","--output",dest="output",default="reconstruct",help="path, where to store the output (default='./reconstruct')")
parser.add_option("","--stdout",dest="stdout",default=None,choices=['circs','lins','reads','multi','test'],help="use to direct chosen type of output (circs, lins, reads, multi) to stdout instead of file")

options,args = parser.parse_args()

# prepare output directory
if not os.path.isdir(options.output):
    os.makedirs(options.output)

# prepare logging system
FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
logging.basicConfig(level=logging.INFO,format=FORMAT,filename=os.path.join(options.output,"run.log"),filemode='w')
logger = logging.getLogger('reconstruct_isoforms.py')
logger.info("reconstruct_isoforms.py {0} invoked as '{1}'".format(__version__," ".join(sys.argv)))


class ExonStorage(object):
    def __init__(self):
        self.sorted_exon_bounds = defaultdict(list)
        self.exons_by_coord = defaultdict( lambda : defaultdict(list) )
        self.logger = logging.getLogger("ExonStorage")

    def load_gtf(self,fname):
        self.logger.info("loading from {0}".format(fname))
        from byo.io.gff import gff_importer
        exon_bounds = defaultdict(set)
        N = 0
        for gff in gff_importer(fname):
            if gff.type != "exon":
                continue
            N += 1
            strand = gff.chrom + gff.sense
            start = gff.start - 1
            end = gff.end
            exon_bounds[strand].add( start )
            exon_bounds[strand].add( end )
            
            self.exons_by_coord[strand][start].append( (start, end) )
            self.exons_by_coord[strand][end].append( (start, end) )
               
        for strand, eb in exon_bounds.items():
            self.sorted_exon_bounds[strand] = sorted(eb)
        
        self.logger.info("done, loaded and sorted {0} exons".format(N))

    def get_intervening_exons(self, chrom, start, end, sense):
        chrom_bounds = self.sorted_exon_bounds[chrom+sense]
        
        start_i = bisect.bisect_left(chrom_bounds, start)
        end_i = bisect.bisect_left(chrom_bounds, end)
        
        iv_bounds = chrom_bounds[start_i:end_i]
        exon_dict = self.exons_by_coord[chrom+sense]
        exons = set()

        for pos in iv_bounds:
            exons |= set(exon_dict[pos])
        
        return sorted(exons)

    
known_exons = ExonStorage()
if options.known_exons:
    known_exons.load_gtf(options.known_exons)

class ReconstructedCircIsoforms(object):
    def __init__(self, multi_event_source):
        self.multi_events = defaultdict(list)

        class MultiEvent(object):
            pass

        def to_coords(column):
            for elem in column.split(','):
                try:
                    start, end = elem.split('-')
                    start = int(start)
                    end = int(end)
                except ValueError:
                    pass
                else:
                    yield (start, end)
            
        for line in multi_event_source:
            if line.startswith("#"):
                continue
            
            parts = line.rstrip().split('\t')
            #print parts
            
            me = MultiEvent()
            me.chrom = parts[0]
            me.start = int(parts[1])
            me.end = int(parts[2])
            me.name = parts[3]
            me.score = parts[4]
            me.sense = parts[5]
            me.read_name = parts[6]

            me.linear = set(to_coords(parts[7]))
            me.unspliced = set(to_coords(parts[9]))

            self.multi_events[me.name].append(me)
        
    def reconstruct(self,circname):
        if not circname in self.multi_events:
            return []
        
        first = self.multi_events[circname][0]
        chrom, start, end, sense = first.chrom, first.start, first.end, first.sense
        candidate_exons = known_exons.get_intervening_exons(chrom, start, end, sense)
        
        print circname, candidate_exons

    def __iter__(self):
        for name in sorted(self.multi_events.keys()):
            if name != "ME:circ_010922":
                continue
            yield self.reconstruct(name)

for rec in ReconstructedCircIsoforms(file(args[0])):
    print rec

