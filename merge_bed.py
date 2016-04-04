#!/usr/bin/env python
import sys,os
from collections import defaultdict
from logging import debug,warning,error,info

from optparse import *

usage = """
%prog 1.bed 2.bed [3.bed] [4.bed] [...] > merged.bed

Merge BED or BED-like files on the genomic coordinates. Deals properly
with find_circ.py output and adds a few extra columns. 
"""

parser = OptionParser(usage=usage)
parser.add_option("-f","--flank",dest="flank",type=int,default=0,help="add flanking nucleotides to define more fuzzy overlap (default=0)")
parser.add_option("-s","--stats",dest="stats",default="",help="write statistics to this file (instead of stderr)")
parser.add_option("-6","--bed6",dest="bed6",default=False,action="store_true",help="ignore all columns except the first six standard BED columns (default=False)")
parser.add_option("-F","--format",dest="format",default="2",choices = ["1","1.2","2"],help="select the find_circ.py outout version (choices=['1','1.2','2'])")
parser.add_option("-V","--verbatim",dest="verbatim",default=False,action="store_true",help="do not attempt to merge all columns. Simply join on coordinates, other columns reported in verbatim")

options,args = parser.parse_args()

from numpy import *

def read_to_hash(fname,ds=0,de=0,flank=0,cover=False):
    #print "loading",fname
    pos = {}
    def src(fname):
        if fname == '-':
            return sys.stdin
        else:
            return file(fname)

    for line in src(fname):
        if line.startswith("#"):
            continue
        line = line.strip()
        parts = line.split('\t')
        if options.bed6:
            parts = parts[:6]

        chrom,start,end,name,score,sense = parts[:6]
        start,end = int(start)+ds,int(end)+de

        #print (chrom,start,end,sense)
        pos[(chrom,start,end,sense)] = parts
        
        if flank:
            for x in xrange(flank):
                pos[(chrom,start-x,end,sense)] = parts
                pos[(chrom,start+x,end,sense)] = parts
                pos[(chrom,start,end-x,sense)] = parts
                pos[(chrom,start,end+x,sense)] = parts
        
        #if cover:
            #for x in xrange
    return pos

N = defaultdict(int)

inputs = [read_to_hash(a,flank=0) for a in args]
names = [os.path.basename(a) for a in args]
shorts = ["in%d" % i for i in range(len(names))]

by_name = dict(zip(shorts,inputs))

merge = {}
support = defaultdict(list)
for name,data in zip(shorts,inputs):
    N[name] = len(data)
    merge.update(data)
    for pos in data:
        support[pos].append(name)

from collections import Counter
comb = Counter([tuple(v) for v in support.values()])

if options.stats:
    sfile = file(options.stats,"w")
else:
    sfile = sys.stderr

for c in sorted(comb.keys()):
    sfile.write("%s\t%d\n" % ("_AND_".join(c),comb[c]))

def consensus_cols(lines,comb):
    samples = []
    counts = defaultdict(int)
    
    def setup_samples(values):
        allsamples = []
        for v in values:
            toadd = v.split(",")
            samples.append(toadd)
            allsamples.extend(toadd)
        return ",".join(sorted(allsamples))
    
    def assign_counts(values):
        for cs,ss in zip(values,samples):
            for samp,count in zip(ss,cs.split(',')):
                counts[samp] += int(count)
            
        res = []
        for k in sorted(counts.keys()):
            res.append(counts[k])
        return ",".join([str(c) for c in res])

    def append_uniq(values):
        v = set()
        for row in values:
            v |= set(row.split(","))

        return ",".join([str(x) for x in sorted(v) if x])

    col_map = { 
        3 : lambda values : ",".join(sorted(values)), # combine identifiers
        4 : lambda values : array(values,dtype=float).sum(), # sum n_reads
        6 : lambda values : array(values,dtype=float).sum(), # sum n_uniq
        7 : lambda values : max([int(x) for x in values]), # max of best_uniq_A
        8 : lambda values : max([int(x) for x in values]), # max of best_uniq_B
        9 : lambda values : array(values,dtype=int).sum(), # sum ov_linear_A
        10 : lambda values : array(values,dtype=int).sum(), # sum ov_linear_B
        11 : setup_samples,
        12 : assign_counts,
        13 : lambda values : min([int(x) for x in values]), # min of edits
        14 : lambda values : min([int(x) for x in values]), # min of anchor_overlap
        15 : lambda values : min([int(x) for x in values]), # min of breakpoints
    }
    
    from itertools import izip_longest
    parts = []
    source = enumerate(izip_longest(*[l for l in lines],fillvalue=""))
    for i,column in source:
        #print i,column
        if i in col_map:
            parts.append(str(col_map[i](column)))
        else:
            parts.append(append_uniq(column))

    return parts

if options.verbatim:
    for pos in merge.keys():
        com = support[pos]
        comstr = "(%s)" % ",".join(com)
        cols = [comstr]
        for name in com:
            cols.append("%s : " % name)
            cols.append("\t".join(by_name[name][pos]))

        print "\t".join(cols)
       
else:
    for pos in merge.keys():
        com = support[pos]
        comstr = "(%s)" % ",".join(com)
        lines = [by_name[name][pos] for name in com]
        cols = [comstr] + consensus_cols(lines,comb)

        print "\t".join(cols)


#for names,inputs in zip(combinations(names


#for circ,line in marv.items():
    #if circ in anna:
        #if len(sys.argv) > 3:
            #print "%s\t%s" % (anna[circ].split('\t')[3],line.split('\t')[3])
        #else:
            #print anna[circ]
        ##print "M",line
        #N['overlap'] += 1        
        #del anna[circ]
    #else:
        #N['input2_not_in_input1'] += 1
    ##print len(anna.keys())
        
#for k,l in anna.items():
    ##if "HEK" in l:
        #print "MISSING\t%s" % l
        #N['input1_not_in_input2'] += 1

#for k in sorted(N.keys()):
    #sys.stderr.write("%s\t%d\n" % (k,N[k]))
        
#found = N['overlap']
#detected = N['unique_input2']
#total = N['unique_input1']
#fp = N['input2_not_in_input1']

#print "#sensitivity %d/%d = %.2f %%" % (found,total,float(found)/total*100)
#print "#FDR %d/%d = %.2f %%" % (fp,detected,float(fp)/detected*100)