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
#logging.basicConfig(level=logging.INFO,format=FORMAT,filename=os.path.join(options.output,"run.log"),filemode='w')
logging.basicConfig(level=logging.DEBUG,format=FORMAT)#,filename=os.path.join(options.output,"run.log"),filemode='w')
logger = logging.getLogger('reconstruct_isoforms.py')
logger.info("reconstruct_isoforms.py {0} invoked as '{1}'".format(__version__," ".join(sys.argv)))


class ExonStorage(object):
    """
    Fast lookup of exons that are inbetween two coordinates
    """
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

class SupportedCircRNA(CircRNA):
    def __init__(self, name, chrom, sense, exon_starts, exon_ends, min_exon_overlap=.85, **kwargs):
        super(CircRNA,self).__init__(name, chrom, sense, exon_starts, exon_ends, (exon_starts[0],exon_starts[0]), **kwargs)
        
        self.min_exon_overlap = min_exon_overlap
        self.exonic_map = {}
        self.junctions = set([ (min(left, right), max(left, right)) for left, right in self.intron_bounds])
        self.reset_counts()
        
        for i,(start,end) in enumerate(self.exon_bounds):
            for x in xrange(start,end):
                self.exonic_map[x] = i

    def reset_counts(self):
        self.exonic_support = defaultdict(int)
        self.junction_support = defaultdict(int)
        
    def splice_support(self,start,end, count=True):
        if not (start,end) in self.junctions:
            return False
        
        if count:
            self.junction_support[(start, end)] += 1
            self.exonic_support[self.exonic_map[start-1]] += 1
            self.exonic_support[self.exonic_map[end]] += 1
        
        return True

    def exon_support(self,start,end, count=True):
        intersect = self.cut(start,end)
        
        if intersect.spliced_length < (end-start) * self.min_exon_overlap:
            return False
        
        if count:
            for x in xrange(start,end):
                if x in self.exonic_map:
                    self.exonic_support[self.exonic_map[x]] += 1
                    break

        return True

    def insert_new_intron(self, left, right):
        new_exon_starts = list(self.exon_starts)
        new_exon_ends = list(self.exon_ends)
        
        start_i = bisect.bisect_left(new_exon_starts, right)
        end_i = bisect.bisect_left(new_exon_ends, left)
        new_exon_starts.insert(start_i, right)
        new_exon_ends.insert(end_i, left)
        
        return SupportedCircRNA("{0}_add_intron_{1}-{2}".format(self.name, left, right), self.chrom, self.sense, new_exon_starts, new_exon_ends)
        
    def adjust_exon_start(self, left, right):
        new_exon_starts = list(self.exon_starts)
        new_exon_ends = list(self.exon_ends)
        
        # find exon by end coordinate
        i = new_exon_ends.index(left)
        new_exon_starts[i] = right
        
        return SupportedCircRNA("{0}_adj_start_{1}:{2}".format(self.name, i, right), self.chrom, self.sense, new_exon_starts, new_exon_ends)

    def adjust_exon_end(self, left, right):
        new_exon_starts = list(self.exon_starts)
        new_exon_ends = list(self.exon_ends)
        
        # find exon by start coordinate
        i = new_exon_starts.index(right)
        new_exon_ends[i] = left
        
        return SupportedCircRNA("{0}_adj_end_{1}:{2}".format(self.name, i, left), self.chrom, self.sense, new_exon_starts, new_exon_ends)

    def skip_exons_between(self, left, right):
        new_exon_starts = list(self.exon_starts)
        new_exon_ends = list(self.exon_ends)

        # find bounding exon indices by coordinates
        i = new_exon_ends.index(left)
        j = new_exon_starts.index(right)
        
        new_exon_ends = new_exon_ends[:i+1] + new_exon_ends[j:]
        new_exon_starts = new_exon_starts[:i+1] + new_exon_starts[j:]
        
        return SupportedCircRNA("{0}_skipped_{1}-{2}".format(self.name, i, j), self.chrom, self.sense, new_exon_starts, new_exon_ends)


    @property
    def support_summary(self):
        print self.junction_support
        print self.exonic_support
        print self.exon_count
        supported_junctions = len(self.junction_support.keys())
        if self.junctions:
            junc_fraction = supported_junctions / float(len(self.junctions))
        else:
            junc_fraction = "n/a"
        
        supported_exons = len(self.exonic_support.keys())
        exon_fraction = supported_exons / float(self.exon_count)

        return "junc_fraction={0} exon_fraction={1}".format(junc_fraction, exon_fraction)



class ReconstructedCircIsoforms(object):
    def __init__(self, multi_event_source):
        self.multi_events = defaultdict(list)
        self.logger = logging.getLogger("ReconstructedCircIsoforms")

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
        
        current_exons_by_bound = defaultdict(list)

        known_bounds = np.array(known_exons.get_intervening_exons(chrom, start, end, sense))
        if len(known_bounds):
            for start, end in known_bounds:
                current_exons_by_bound[start].append( (start, end) )
                current_exons_by_bound[end].append( (start, end) )

            starts, ends = known_bounds.transpose()
            isoforms = [SupportedCircRNA("{0}_known".format(circname),chrom, sense, starts, ends)]
            self.logger.info("starting reconstruction of {circname} from {n} known exons".format(circname=circname, n=len(known_bounds) ))
        else:
            # no known exons in this region? 
            # Assume single exon circRNA
            starts = [first.start]
            ends = [first.end]
            current_exons_by_bound[first.start].append( (first.start, first.end) )
            current_exons_by_bound[first.end].append( (first.start, first.end) )
            
            isoforms = [SupportedCircRNA("{0}_denovo".format(circname),chrom, sense, starts, ends)]
            self.logger.info("starting reconstruction of {circname} de novo".format(circname=circname) )


        # TODO: 
        # compute a square support matrix (including junctions and coverage) 
        # between reads and isoforms to:
        # a) identify the best matching isoform to start with
        # b) prune redundant isoforms in the end
        for me in self.multi_events[circname]:
            for left, right in me.linear:
                supported = np.array([I.splice_support(left, right) for I in isoforms])
                if not supported.any():
                    best = isoforms[-1] # TODO: sort by best support, including exon overlap
                    
                    has_left = left in best.exon_ends
                    has_right = right in best.exon_starts
                    
                    if not has_left and not has_right:
                        self.logger.info("discovered new intron {left}-{right}".format(**locals()))
                        new_iso = best.insert_new_intron(left, right)
                        assert new_iso.splice_support(left, right)
                        isoforms.append(new_iso)
                        
                    elif not has_right:
                        best = isoforms[-1]
                        self.logger.info("discovered alternative exon start {left}".format(left=left))
                        new_iso = best.adjust_exon_start(left, right)
                        assert new_iso.splice_support(left, right)
                        isoforms.append(new_iso)

                    elif not has_left:
                        best = isoforms[-1]
                        self.logger.info("discovered alternative exon end {left}".format(right=right))
                        new_iso = best.adjust_exon_end(left, right)
                        assert new_iso.splice_support(left, right)
                        isoforms.append(new_iso)

                    else:
                        best = isoform[-1]
                        self.logger.info("discovered skipped exon {left}-{right}".format(left=left, right=right))
                        new_iso = best.skip_exons_between(left, right)
                        assert new_iso.splice_support(left, right)
                        isoforms.append(new_iso)
                        

            for start, end in me.unspliced:
                supported = np.array([I.exon_support(start,end) for I in isoforms])
                if not supported.any():
                    print "Need new exon or retained intron"
                    
        return isoforms
 


    def __iter__(self):
        for name in sorted(self.multi_events.keys()):
            #if name != "ME:circ_010922":
                #continue
            print name
            yield self.reconstruct(name)

for isoform_set in ReconstructedCircIsoforms(file(args[0])):
    for iso in isoform_set:
        print iso.support_summary

