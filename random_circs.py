#!/usr/bin/env python
from byo.gene_model import transcripts_from_GTF, Transcript
from collections import defaultdict
import optparse
import logging
import sys
import numpy as np

FORMAT = '%(asctime)-20s\t%(levelname)s\t%(name)s\t%(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
log = logging.getLogger("random_circs.py")

usage = """
   cat exon_nums | %prog [options] > simulated_circRNA_models.ucsc
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
parser.add_option("","--seed",dest="seed",type=int,default=0,help="seed for pseudo random number generator (default=random)")
options,args = parser.parse_args()

tx_by_exon_count = defaultdict(list)
exons_used = set()
n_max = 0

log.info("sorting known transcripts")
for tx in transcripts_from_GTF(args[0]):
    tx_by_exon_count[tx.exon_count].append(tx)
    n_max = max(n_max, tx.exon_count)
    
if options.seed:
    log.info("seeding random number generator with {0}".format(options.seed))
    np.random.seed(options.seed)

n_sim = 0
log.info("generating random circ structures")

def rand_start(max_n):
    if max_n == 1:
        return 1
    return np.random.randint(1, max_n)
    
for line in sys.stdin:
    n_exons = int(line)
    
    candidates = []
    for n in range(n_exons+2,n_max):
        candidates.extend(tx_by_exon_count[n])
    
    log.debug("{0} candidate transcripts for {1} exon circ".format(len(candidates), n_exons))
    candidate_i = list(np.random.permutation(len(candidates)))
    
    circ = None
    while not circ and candidate_i:
        cand_i = candidate_i.pop()
        cand = candidates[cand_i]

        n_first = rand_start(cand.exon_count-n_exons-1)
        n_last = n_first + n_exons
        
        exon_starts = cand.exon_starts[n_first:n_last]
        exon_ends = cand.exon_ends[n_first:n_last]
        
        # test if exons have already been used?
        exon_keys = [(cand.chrom, cand.sense, s, e) for s,e in zip(exon_starts, exon_ends)]
        
        is_used = False
        for key in exon_keys:
            if key in exons_used:
                is_used = True
                break

        if is_used:
            continue

        # remember the used exons
        exons_used |= set(exon_keys)
        n_sim += 1
        name = "syncirc{0}".format(n_sim)
        circ = Transcript(name, cand.chrom, cand.sense, exon_starts, exon_ends, (exon_starts[0], exon_starts[0]), system = options.system)
        
    if not circ:
        log.warning("exhausted all candidates for {0}!".format(n_exons))
        continue
    else:
        print circ
