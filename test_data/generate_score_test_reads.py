#!/usr/bin/env python

from byo.systems import hg19

RG = hg19.get_refGenes()
elavl1 = RG['ENST00000407627.2']
mbnl1 = RG['ENST00000498502.1']

elavl1_exons = list(elavl1.exons)
mbnl1_exons = list(mbnl1.exons)


ex2, ex3, ex4 = [hg19.genome.get_oriented(exon.chrom, exon.start, exon.end, exon.sense) for exon in elavl1_exons[2:5]]
mx2, mx3, mx4 = [hg19.genome.get_oriented(exon.chrom, exon.start, exon.end, exon.sense) for exon in mbnl1_exons[2:5]]
print ">pos_ctrl_100\n{0}\n".format(ex2[-100:]+ex3 + ex2[:100])
print ">pos_ctrl_20\n{0}\n".format(ex2[-100:]+ex3 + ex2[:20])
print ">pos_ctrl_18\n{0}\n".format(ex2[-100:]+ex3 + ex2[:18])
print ">pos_ctrl_15\n{0}\n".format(ex2[-100:]+ex3 + ex2[:15])

print ">pos_ctrl_for_fusion\n{0}\n".format(ex3[-100:]+ex2)
print ">neg_ctrl_fusion_100\n{0}\n".format(ex3[-100:]+ex2+mx3[:100])
print ">neg_ctrl_trans_100\n{0}\n".format(ex3[-100:]+ex2+ex3+ex4[:100])
print ">neg_ctrl_trans_20\n{0}\n".format(ex3[-100:]+ex2+ex3+ex4[:20])
print ">neg_ctrl_trans_18\n{0}\n".format(ex3[-100:]+ex2+ex3+ex4[:18])
print ">neg_ctrl_trans_15\n{0}\n".format(ex3[-100:]+ex2+ex3+ex4[:15])





