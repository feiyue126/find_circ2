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
        super(SupportedCircRNA,self).__init__(name, chrom, sense, exon_starts, exon_ends, (exon_starts[0],exon_starts[0]), **kwargs)
        
        self.logger = logging.getLogger("SupportedCircRNA")
        self.min_exon_overlap = min_exon_overlap
        self.exonic_map = {}
        self.junctions = set([ (min(left, right), max(left, right)) for left, right in self.intron_bounds])
        self.exon_starts_set = set(self.exon_starts) 
        self.exon_ends_set = set(self.exon_ends)
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

    def compatibility_score(self, me):
        
        max_junc_score = 12. * len(me.linear) + 2.
        junc_score = \
            10. * len( me.linear & self.junctions) + \
            len(me.exon_starts_set & self.exon_starts_set) + \
            len(me.exon_ends_set & self.exon_ends_set)

        junc_ratio = junc_score/max_junc_score
        
        max_exon_score = 0.
        exon_score = 0
        for start,end in me.unspliced:
            L = end - start
            max_exon_score += L * self.min_exon_overlap
            
            overlap = 0
            for x in xrange(start,end):
                if x in self.exonic_map:
                    overlap += 1

            # allow a little bit of misaligned coverage (default=15%), so max-score is 0.85 * L
            exon_score += min(self.min_exon_overlap*L, overlap )
            
        #self.logger.debug('junc_ratio={0} max_exon_score={1} exon_score={2}'.format(junc_ratio, max_exon_score, exon_score))
        
        if max_exon_score:
            return 0.75 * junc_ratio + 0.25 * exon_score/max_exon_score
        else:
            return junc_ratio
        
        
    @property
    def support_summary(self):
        #print self.junction_support
        #print self.exonic_support
        #print self.exon_count
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
            me.exon_starts = [me.start,] + [l[1] for l in me.linear]
            me.exon_ends = [l[0] for l in me.linear] + [me.end,]
            me.exon_starts_set = set(me.exon_starts)
            me.exon_ends_set = set(me.exon_ends)

            me.exon_bound_lookup = {}
            for start, end in zip(me.exon_starts, me.exon_ends):
                me.exon_bound_lookup[start] = end
                me.exon_bound_lookup[end] = start
            
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
            isoforms = [SupportedCircRNA("{0}:known".format(circname),chrom, sense, starts, ends)]
            self.logger.info("starting reconstruction of {circname} from {n} known exons".format(circname=circname, n=len(known_bounds) ))
        else:
            # no known exons in this region? 
            # Assume single exon circRNA
            starts = [first.start]
            ends = [first.end]
            current_exons_by_bound[first.start].append( (first.start, first.end) )
            current_exons_by_bound[first.end].append( (first.start, first.end) )
            
            isoforms = [SupportedCircRNA("{0}:denovo".format(circname),chrom, sense, starts, ends)]
            self.logger.info("starting reconstruction of {circname} de novo".format(circname=circname) )

        all_multi_events = self.multi_events[circname]
        self.logger.debug("processing {0} multi events".format(len(all_multi_events)) )
        
        # compute a square support matrix (including junctions and coverage) 
        # between reads and isoforms to:
        # a) identify the best matching isoform to start with
        # b) prune redundant isoforms in the end

        matrix = np.array([[I.compatibility_score(me) for me in all_multi_events] for I in isoforms])
        #best support (by any isoform) for each me is matrix.max(axis=0)
        to_process = list((matrix.max(axis=0) < 1.).nonzero()[0])
        #print "the following me's are not fully compatible with any isoform in the list", to_process
        
        while len(to_process):
            # pick the first incompatible me
            i = to_process.pop(0)
            me = all_multi_events[i]
            
            # and find the best matching current isoform to start with
            ties = (matrix[:,i] == matrix[:,i].max() ).nonzero()[0]
            logger.debug("selecting candidates from {0}".format(matrix[:,i]) )
            
            # in case of ties, pick the one with more junctions to avoid fragmentation
            scores = [(isoforms[ind].exon_count, ind) for ind in ties]
            logger.debug("ties={0}".format(ties))
            logger.debug("scores={0}".format(scores))
            j = sorted(scores, reverse=True)[0][1]
            
            best = isoforms[j]
            logger.debug("best matched isoform {0}".format(best) )
            
            self.logger.debug("selected incompatible me {0} and best-matched isoform {1} at compatibility_score {2}".format(i,j, matrix[j][i]) )
            
            # first take care of completely unknown junctions
            for left, right in (me.linear - best.junctions):
                if (left in best.exon_ends_set) and (right in best.exon_starts_set):
                    self.logger.info("  discovered skipped exon {left}-{right}".format(left=left, right=right))
                    new_iso = best.skip_exons_between(left, right)
                    assert new_iso.splice_support(left, right)
                else:
                    self.logger.info("  discovered new intron {left}-{right}".format(**locals()))
                    new_iso = best.insert_new_intron(left, right)
                    assert new_iso.splice_support(left, right)

                best = new_iso

            # next, take care of alternative 5' and 3' end positions.
            for start in (me.exon_starts_set - best.exon_starts_set):
                end = me.exon_bound_lookup[start]
                self.logger.info("  discovered alternative exon start {left}".format(left=left))
                new_iso = best.adjust_exon_start(start, end)
                assert new_iso.splice_support(start, end)
                best = new_iso

            for end in (me.exon_ends_set - best.exon_ends_set):
                start = me.exon_bound_lookup[end]
                self.logger.info("  discovered alternative exon end {left}".format(right=right))
                new_iso = best.adjust_exon_end(start, end)
                assert new_iso.splice_support(start, end)
                best = new_iso

            # finally, if there is unspliced coverage in an intronic region, remove this intron
            for start, end in me.unspliced:
                if not best.exon_support(start,end):
                    #print "Need new exon or retained intron"
                    self.logger.info("  removing intron with unspliced coverage from {start}-{end}".format(start=start, end=end) )
                    best = best.remove_intron_overlapping(start, end)

            new_compat = best.compatibility_score(me)
            if new_compat < 1.:
                self.logger.warning("  unsatisfyable multievent. Ignoring {0}".format(me))
                all_multi_events.remove(me)
            else:
                isoforms.append(best)
                
            # recompute matrix
            matrix = np.array([[I.compatibility_score(me) for me in all_multi_events] for I in isoforms])
            to_process = list((matrix.max(axis=0) < 1.).nonzero()[0])
            self.logger.debug(" recomputed compatibility matrix with {0} incompatible ME's left".format(len(to_process) ) )
        
        # prune non-essential isoforms and sort by ME support
        me_support = (matrix == 1).sum(axis=1)
        me_matched= np.zeros(len(all_multi_events), dtype=bool)
        logger.debug("me_support {0}".format(me_support) )
        minimal_set = []
        for i in me_support.argsort()[::-1]:
            matches = (matrix[i] == 1.)
            # are more me's matched if we accept this isoform?
            if (me_matched | matches).sum() > me_matched.sum():
                # record number of supporting/compatible ME's as score
                isoforms[i].score = me_support[i]
                minimal_set.append(isoforms[i])
                me_matched |= matches
            
        self.logger.debug("---> returning minimal set of {0} isoforms compatible with all {1} multi-events".format(len(minimal_set), len(all_multi_events)) )
        return minimal_set
 

    def __iter__(self):
        for name in sorted(self.multi_events.keys()):
            #if name != "ME:circ_010922":
                #continue
            #print name
            yield self.reconstruct(name)

for isoform_set in ReconstructedCircIsoforms(file(args[0])):
    for iso in isoform_set:
        print iso
