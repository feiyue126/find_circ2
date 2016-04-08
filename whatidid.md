# new BWA MEM support tested on known spliced reads from HEK293 find_circ.py (bowtie2) run:
```
    cat hek_test_out/spliced_reads.fa | sed 's/ /_/g' > test_bwa_mem_reads.fa
    grep CIRCULAR hek_test_out/splice_sites.bed > bt2_circs.bed
    bwa mem -t16 -k 15 -T 1 /data/rajewsky/indices/hg19_bwa_0.7.12-r1039/hg19.fa test_bwa_mem_reads.fa > bwa_mem_test.sam

    cat bwa_mem_test.sam | ../find_circ.py -S hg19 --bwa-mem > bwamem.spliced.bed
    grep CIRCULAR bwamem.spliced.bed > bm_circs.bed
    ../cmp_bed.py bt2_circs.bed bm_circs.bed 
```


# making test 150nt reads that should align in three segments
```
    grep ENST00000407627 ~/pcp/systems/hg19/annotation/wgEncodeGencodeBasicV17.ucsc | python simulate_reads.py > elavl_test_reads.fa
    bwa mem -t16 -k 15 -T 1 /data/rajewsky/indices/hg19_bwa_0.7.12-r1039/hg19.fa elavl_test_reads.fa > bwa_mem_triplet.sam
    cat bwa_mem_triplet.sam | ../find_circ.py -S hg19 --bwa-mem --debug
```
