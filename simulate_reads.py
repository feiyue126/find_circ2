#!/usr/bin/env python
from byo.gene_model import transcripts_from_UCSC, CircRNA
from byo import rev_comp
import os,sys
import optparse
import numpy as np

usage = """
   cat circ_models.ucsc | %prog [options] > interlaced_paired_end_reads.fa
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-S","--system",dest="system",type=str,default="",help="model system database (optional! Requires byo library.)")
parser.add_option("","--n-frags",dest="n_frags",type=int,default=100,help="number of fragments to simulate (default=100)")
parser.add_option("","--frag-len",dest="frag_len",type=int,default=350,help="fragment length to simulate (default=350)")
parser.add_option("","--read-len",dest="read_len",type=int,default=100,help="read length to simulate (default=100)")
options,args = parser.parse_args()

import importlib
system = importlib.import_module("byo.systems.{options.system}".format(**locals()))

def test_str(mate):
    ori = mate.origin
    #print "mate origin",mate.origin, mate.origin_spliced, mate.map_to_exon(mate.origin_spliced+1)
    seg_bounds = np.roll(mate.exon_bounds, -mate.map_to_exon(mate.origin_spliced+1), axis=0)
    #print mate.name, seg_bounds

    last_s, last_e = seg_bounds[0]

    parts = []
    for i,(s,e) in enumerate(seg_bounds[::mate.dir]):
        if i == 0:
            # first segment: yield origin
            parts.append("O:{chrom}:{ori}:{sense};M:{m}".format(chrom=mate.chrom, ori=mate.origin, sense=mate.sense, m=e-s) )
        elif s > last_s:
            # we move forward: linear splicing
            parts.append("LS:{start}:{end};M:{m}".format(start=last_e - ori, end=s - ori, m=e-s) )
        elif s < last_s:
            # we move backward: circRNA backsplicing
            parts.append("CS:{start}:{end};M:{m}".format(start=s - ori, end=last_e - ori, m=e-s) )
    
        last_s, last_e = s,e

    return ";".join(parts)

for circ in transcripts_from_UCSC(sys.stdin, system=system, tx_type=CircRNA):

    L = circ.spliced_length
    rnd_pos = np.random.randint(0, L, options.n_frags)
    #rnd_pos = [100,600,900]
    for i,frag_start in enumerate(rnd_pos):
        
        m1_start, m1_end = frag_start, frag_start + options.read_len
        m2_start, m2_end = frag_start + options.frag_len - options.read_len, frag_start + options.frag_len
        
        m1_g_start, m1_g_end = circ.map_from_spliced(m1_start), circ.map_from_spliced(m1_end)
        m2_g_start, m2_g_end = circ.map_from_spliced(m2_start), circ.map_from_spliced(m2_end)
        
        mate1 = circ.cut(m1_g_start, m1_g_end)
        mate2 = circ.cut(m2_g_start, m2_g_end)
        
        #print "start_pos",frag_start
        #print "circ"
        #print circ
        #print "rel. coordinates"
        #print L
        #print "m1:",m1_start,m1_end
        #print "m2:",m2_start,m2_end
        
        #print "gen. coordinates"
        #print "m1:",m1_g_start,m1_g_end
        #print "m2:",m2_g_start,m2_g_end
        
        #print "resulting chains"
        #print mate1
        #print mate2
        #print mate1.spliced_length
        #print mate2.spliced_length
        
        assert mate1.spliced_length == options.read_len
        assert mate2.spliced_length == options.read_len
        
        # simulate forward/reverse mate pairs
        mate1_seq = mate1.spliced_sequence
        mate2_seq = rev_comp(mate2.spliced_sequence)
        
        m1_str = test_str(mate1)
        m2_str = test_str(mate2)

        read_name = "sim_{circ.name}_{i}___{m1_str}|{m2_str}".format(**locals())
        print ">{read_name}\n{mate1_seq}\n>{read_name}\n{mate2_seq}".format(**locals())
