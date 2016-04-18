#!/usr/bin/env python

__version__ = "0.90"
__author__ = "Marvin Jens"
__credits__ = ["Marvin Jens","Marcel Schilling","Petar Glazar","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marvin.jens@mdc-berlin.de"

from byo.gene_model import transcripts_from_UCSC, CircRNA, Transcript
#from byo import rev_comp
import os,sys
import optparse
import numpy as np
import bisect
import logging
from collections import defaultdict
#from gzip import GzipFile

usage = """
   %prog [options] multi_events.tsv
"""

parser = optparse.OptionParser(usage=usage)
#parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
#parser.add_option("-G","--genome",dest="genome",type=str,default="",help="path to genome FASTA file")
parser.add_option("","--known-exons",dest="known_exons",default="",help="GTF file with known exons")
parser.add_option("","--max-isoforms",dest="max_isoforms",type=int,default=10,help="maximal number of candidate isoforms to investigate (default=10)")
parser.add_option("","--max-size",dest="max_size",type=int,default=100,help="maximal size of genomic region to investigate in kb, to prevent excessive RAM usage (default=100kb)")
parser.add_option("","--max-exon-size",dest="max_exon_size",type=int,default=10000,help="maximal size of an exon in nt to prevent excessive RAM usage (default=10000)")
#parser.add_option("","--n-frags",dest="n_frags",type=int,default=100,help="number of fragments to simulate (default=100)")
parser.add_option("","--self-test",dest="self_test",action="store_true",default=False, help="if set, perform unit-tests instead of running on input data")
parser.add_option("","--debug",dest="debug",default=False,action="store_true",help="Activate LOTS of debug output")
#parser.add_option("","--frag-len",dest="frag_len",type=int,default=350,help="fragment length to simulate (default=350)")
#parser.add_option("","--read-len",dest="read_len",type=int,default=100,help="read length to simulate (default=100)")
parser.add_option("-o","--output",dest="output",default="reconstruct",help="path, where to store the output (default='./reconstruct')")
parser.add_option("","--stdout",dest="stdout",default=None,choices=['circs','lins','reads','multi','test'],help="use to direct chosen type of output (circs, lins, reads, multi) to stdout instead of file")

options,args = parser.parse_args()

# prepare output directory
if not os.path.isdir(options.output):
    os.makedirs(options.output)

# prepare logging system
if options.debug:
    lvl = logging.DEBUG
else:
    lvl = logging.INFO
FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
logging.basicConfig(level=lvl,format=FORMAT,filename=os.path.join(options.output,"reconstruct.log"),filemode='w')
logger = logging.getLogger('reconstruct_isoforms.py')
logger.info("reconstruct_isoforms.py {0} invoked as '{1}'".format(__version__," ".join(sys.argv)))

# prepare output files
iso_file = file(os.path.join(options.output,"reconstructed.ucsc"),"w")
psi_file  = file(os.path.join(options.output,"psi.bed"),"w")

# redirect output to stdout, if requested
if options.stdout:
    varname = "{0}_file".format(options.stdout)
    out_file = globals()[varname]
    out_file.write('# redirected to stdout\n')
    logger.info('redirected {0} to stdout'.format(options.stdout))
    globals()[varname] = sys.stdout


class ExonStorage(object):
    """
    Fast lookup of exons that are inbetween two coordinates
    """
    def __init__(self):
        self.sorted_exon_bounds = defaultdict(list)
        self.exons_by_coord = defaultdict( lambda : defaultdict(set) )
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
            
            #self.exons_by_coord[strand][start].add( (start, end) )
            #self.exons_by_coord[strand][end].add( (start, end) )
               
        for strand, eb in exon_bounds.items():
            self.sorted_exon_bounds[strand] = sorted(eb)
            
        #for strand in self.exons_by_coord.keys():
            #for coord, exon_list in self.exons_by_coord[strand].keys():
                #self.exons_by_coord[strand][coord] = sorted(set(exon_list))
        
        self.logger.info("done, loaded and sorted {0} exons".format(N))


    def add(self, chrom, start, end, sense):
        strand = chrom+sense
        eb = (start, end)
        
        #sorted_coords = self.exons_by_coord[strand]
        #sorted_coords[start].add( eb )
        #sorted_coords[end].add( eb )
        
        sorted_bounds = self.sorted_exon_bounds[strand]
        sorted_bounds.insert(bisect.bisect_left(sorted_bounds, eb ), eb )
        

    def get_intervening_exons(self, chrom, start, end, sense):
        chrom_bounds = self.sorted_exon_bounds[chrom+sense]
        
        start_i = bisect.bisect_left(chrom_bounds, (start, start) )
        end_i = bisect.bisect_right(chrom_bounds, (end, end) )
        
        iv_bounds = chrom_bounds[start_i:end_i]
        return iv_bounds
    
known_exons = ExonStorage()
if options.known_exons:
    known_exons.load_gtf(options.known_exons)

class SupportedCircRNA(CircRNA):
    def __init__(self, name, chrom, sense, exon_starts, exon_ends, min_exon_overlap=.75, **kwargs):
        super(SupportedCircRNA,self).__init__(name, chrom, sense, exon_starts, exon_ends, (exon_starts[0],exon_starts[0]), **kwargs)
        
        self.logger = logging.getLogger("SupportedCircRNA")
        self.min_exon_overlap = min_exon_overlap
        self.exonic_map = {}
        self.junctions = set([ (min(left, right), max(left, right)) for left, right in self.intron_bounds])
        self.exon_starts_set = set(self.exon_starts) 
        self.exon_ends_set = set(self.exon_ends)
        self.reset_counts()

        for i,(start,end) in enumerate(self.exon_bounds):
            if (end - start) > options.max_exon_size:
                self.logger.warning("exon{0} of circRNA {1} is {2} kb, exceeding --max-exon-size. Disabling exonic_map".format(i, name, (end-start)/1000.) )
            else:
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
    
    def remove_intron_overlapping(self, left, right):
        #self.logger.debug("remove_intron_overlapping(): initial exon_starts={0} exon_ends={1}".format(self.exon_starts, self.exon_ends) )
        start_i = max(bisect.bisect_left(self.exon_starts, right),1)
        end_i = min(bisect.bisect_left(self.exon_ends, left)+1,self.exon_count-1)
        #self.logger.debug("remove_intron_overlapping(): start_i={0} end_i={1}".format(start_i, end_i) )
        
        drop_left= self.exon_ends[start_i-1]
        drop_right= self.exon_starts[end_i]
        #self.logger.debug("remove_intron_overlapping(): drop_left={0} drop_right={1}".format(drop_left, drop_right) )
        
        new_exon_starts = list(self.exon_starts[:start_i]) + list(self.exon_starts[end_i+1:])
        new_exon_ends = list(self.exon_ends[:start_i-1]) + list(self.exon_ends[end_i:])
        
        #self.logger.debug("remove_intron_overlapping(): new exon_starts={0} exon_ends={1}".format(new_exon_starts, new_exon_ends) )
        return SupportedCircRNA("{0}_drop_intron_{1}-{2}".format(self.name, drop_left, drop_right), self.chrom, self.sense, new_exon_starts, new_exon_ends)
        
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
        
        return SupportedCircRNA("{0}_skipped_exons_{1}-{2}".format(self.name, i+1, j-1), self.chrom, self.sense, new_exon_starts, new_exon_ends)

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

        if options.debug:
            self.logger.debug('junc_ratio={0} max_exon_score={1} exon_score={2}'.format(junc_ratio, max_exon_score, exon_score))
        
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

class MultiEvent(object):
    def __init__(self, name, chrom, start, end, score, sense, read_name = "readname", linear = [], unspliced = []):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.score = score
        self.sense = sense
        self.read_name = read_name
        self.linear = set(linear)
        self.unspliced = set(unspliced)
        
        self.exon_starts = [self.start,] + [l[1] for l in self.linear]
        self.exon_ends = [l[0] for l in self.linear] + [self.end,]
        self.exon_starts_set = set(self.exon_starts)
        self.exon_ends_set = set(self.exon_ends)

        self.exon_bound_lookup = {}
        for start, end in zip(self.exon_starts, self.exon_ends):
            self.exon_bound_lookup[start] = end
            self.exon_bound_lookup[end] = start

    def __str__(self):
        return "MultiEvent({self.chrom}, start={self.start}, end={self.end}, score={self.score}, sense={self.sense}, read_name={self.read_name}, linear={self.linear}, unspliced={self.unspliced})".format(self=self)

class ReconstructedCircIsoforms(object):
    def __init__(self, multi_event_source, known_exon_storage = known_exons):
        self.multi_events = defaultdict(list)
        self.logger = logging.getLogger("ReconstructedCircIsoforms")
        self.known_exon_storage = known_exon_storage
        for me in multi_event_source:
            self.multi_events[me.name].append(me)

        # TODO: percent spliced in metric output
        self.psi = {}

    def reconstruct(self,circname):
        if not circname in self.multi_events:
            return []
        
        first = self.multi_events[circname][0]
        chrom, start, end, sense = first.chrom, first.start, first.end, first.sense
        
        current_exons_by_bound = defaultdict(list)

        known_bounds = np.array(self.known_exon_storage.get_intervening_exons(chrom, start, end, sense))
        if len(known_bounds):
            for start, end in known_bounds:
                current_exons_by_bound[start].append( (start, end) )
                current_exons_by_bound[end].append( (start, end) )

            starts, ends = known_bounds.transpose()
            initial = SupportedCircRNA("{0}:known".format(circname),chrom, sense, starts, ends)
            self.logger.info("starting reconstruction of {circname} from {n} known exons".format(circname=circname, n=len(known_bounds) ))
        else:
            # no known exons in this region? 
            # Assume single exon circRNA
            starts = [first.start]
            ends = [first.end]
            current_exons_by_bound[first.start].append( (first.start, first.end) )
            current_exons_by_bound[first.end].append( (first.start, first.end) )
            
            initial = SupportedCircRNA("{0}:denovo".format(circname),chrom, sense, starts, ends)
            self.logger.info("starting reconstruction of {circname} de novo".format(circname=circname) )

        if (initial.end - initial.start)/1000. > options.max_size:
            self.logger.error("circRNA {0} exceeds --max-size in genomic size".format(circname))
            return []
        
        isoforms = [initial]
        all_multi_events = self.multi_events[circname]
        self.logger.debug("processing {0} multi events".format(len(all_multi_events)) )
        self.logger.debug("initial isoform: '{0}'".format(isoforms[0]) )
        # compute a square support matrix (including junctions and coverage) 
        # between reads and isoforms to:
        # a) identify the best matching isoform to start with
        # b) prune redundant isoforms in the end

        matrix = np.array([[I.compatibility_score(me) for me in all_multi_events] for I in isoforms])
        #best support (by any isoform) for each me is matrix.max(axis=0)
        to_process = list((matrix.max(axis=0) < 1.).nonzero()[0])
        #print "the following me's are not fully compatible with any isoform in the list", to_process
        
        while len(to_process):
            if len(isoforms) > options.max_isoforms:
                self.logger.error("exhausted maximal number of candidate isoforms on circRNA {0}".format(circname))
                return []
            # pick the first incompatible me
            i = to_process.pop(0)
            me = all_multi_events[i]
            
            # and find the best matching current isoform to start with
            ties = (matrix[:,i] == matrix[:,i].max() ).nonzero()[0]
            self.logger.debug("selecting candidates from {0}".format(matrix[:,i]) )
            
            # in case of ties, pick the one with more junctions to avoid fragmentation
            scores = [(isoforms[ind].exon_count, ind) for ind in ties]
            if options.debug:
                self.logger.debug("ties={0}".format(ties))
                self.logger.debug("scores={0}".format(scores))
            j = sorted(scores, reverse=True)[0][1]
            
            best = isoforms[j]
            if options.debug:
                self.logger.debug("best matched isoform {0}".format(best) )
                self.logger.debug("selected incompatible me {0} (readname={3}) and best-matched isoform {1} ({4}) at compatibility_score {2}".format(i,j, matrix[j][i], me.read_name, best.name) )
            
            # first take care of completely unknown junctions
            for left, right in (me.linear - best.junctions):
                if (left in best.exon_ends_set) and (right in best.exon_starts_set):
                    self.logger.info("  discovered skipped exon {left}-{right}".format(left=left, right=right))
                    new_iso = best.skip_exons_between(left, right)
                    assert new_iso.splice_support(left, right, count=False)
                else:
                    self.logger.info("  discovered new intron {left}-{right}".format(**locals()))
                    new_iso = best.insert_new_intron(left, right)
                    assert new_iso.splice_support(left, right, count=False)

                best = new_iso

            # next, take care of alternative 5' and 3' end positions.
            for start in (me.exon_starts_set - best.exon_starts_set):
                end = me.exon_bound_lookup[start]
                self.logger.info("  discovered alternative exon start {left}".format(left=left))
                new_iso = best.adjust_exon_start(start, end, count=False)
                assert new_iso.splice_support(start, end, count=False)
                best = new_iso

            for end in (me.exon_ends_set - best.exon_ends_set):
                start = me.exon_bound_lookup[end]
                self.logger.info("  discovered alternative exon end {left}".format(right=right))
                new_iso = best.adjust_exon_end(start, end, count=False)
                assert new_iso.splice_support(start, end, count=False)
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

def multi_events_from_file(fname):
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
        
    for line in file(fname):
        if line.startswith("#"):
            continue
        
        parts = line.rstrip().split('\t')
        #print parts
        
        chrom, start, end, name, score, sense = parts[:6]
        me = MultiEvent(
            name, chrom, int(start), int(end), float(score), sense, 
            read_name=parts[6], 
            linear = set(to_coords(parts[7])), 
            unspliced = set(to_coords(parts[9]))
        )
        yield me
    

def test_reconstruction():
    print "running self-test 'double_exon'"
    multi_events = [
        MultiEvent("test_double_exon", "chrNA", 10, 100, 1, '+', read_name='test1_perfect_cov_exon1', linear = [(30,70)], unspliced = [(11,29)]),
        MultiEvent("test_double_exon", "chrNA", 10, 100, 1, '+', read_name='test1_perfect_cov_exon2', linear = [(30,70)], unspliced = [(71,99)]),
        MultiEvent("test_double_exon", "chrNA", 10, 100, 1, '+', read_name='test1_bound_cov_exon1', linear = [(30,70)], unspliced = [(10,30)]),
        MultiEvent("test_double_exon", "chrNA", 10, 100, 1, '+', read_name='test1_bound_cov_exon2', linear = [(30,70)], unspliced = [(70,100)]),
        MultiEvent("test_double_exon", "chrNA", 10, 100, 1, '+', read_name='test1_exceed_cov_exon1', linear = [(30,70)], unspliced = [(8,32)]),
        MultiEvent("test_double_exon", "chrNA", 10, 100, 1, '+', read_name='test1_exceed_cov_exon2', linear = [(30,70)], unspliced = [(65,102)]),
    ]
    for isoform_set in ReconstructedCircIsoforms(multi_events):
        for iso in isoform_set:
            print iso
    
    print "running self-test 'known_triple_exon_intron_retention'"
    test_exons = ExonStorage()
    test_exons.add('chrNA', 10, 30, '+')
    test_exons.add('chrNA', 70, 100, '+')
    test_exons.add('chrNA', 150, 200, '+')

    multi_events = [
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_perfect_cov_exon1', linear = [(30,70)], unspliced = [(11,29)]),
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_perfect_cov_exon2', linear = [(30,70)], unspliced = [(71,99)]),
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_perfect_cov_exon2', linear = [(100,150)], unspliced = [(71,99)]),
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_IR_intron2', linear = [(30,70)], unspliced = [(90,120)]),
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_IR_intron1', linear = [(100,150)], unspliced = [(30,70)]),
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_IR_introns1_2', linear = [], unspliced = [(35,130)]),
        MultiEvent("test_known_triple_ir", "chrNA", 10, 200, 1, '+', read_name='test2_IR', linear = [], unspliced = [(40,90)]),
    ]
    for isoform_set in ReconstructedCircIsoforms(multi_events, known_exon_storage=test_exons):
        for iso in isoform_set:
            print iso

    print "running self-test 'skipped_exon'"
    test_exons = ExonStorage()
    test_exons.add('chrNA', 10, 30, '+')
    test_exons.add('chrNA', 70, 100, '+')
    test_exons.add('chrNA', 150, 200, '+')

    multi_events = [
        MultiEvent("test_known_triple_AS", "chrNA", 10, 200, 1, '+', read_name='test2_perfect_cov_exon1', linear = [(30,70)], unspliced = [(11,29)]),
        MultiEvent("test_known_triple_AS", "chrNA", 10, 200, 1, '+', read_name='test2_perfect_cov_exon2', linear = [(30,70)], unspliced = [(71,99)]),
        MultiEvent("test_known_triple_AS", "chrNA", 10, 200, 1, '+', read_name='test2_perfect_complete_splice', linear = [(30,70),(100,150)], unspliced = [(71,99)]),
        MultiEvent("test_known_triple_AS", "chrNA", 10, 200, 1, '+', read_name='test2_exon2_skipped', linear = [(30,150)], unspliced = [(11,29)]),
    ]
    for isoform_set in ReconstructedCircIsoforms(multi_events, known_exon_storage=test_exons):
        for iso in isoform_set:
            print iso


    
if options.self_test:
    test_reconstruction()
    
else:
    for isoform_set in ReconstructedCircIsoforms(multi_events_from_file(args[0])):
        for iso in isoform_set:
            iso_file.write(str(iso) + '\n')
