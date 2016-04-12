## Testing the --test feature

# Using Marcel's C.elegans reads

```
    # in find_circ2/test_data
    PREFIX=/data/BIO3/home/mschilli/repo/global/external/marvin/find_circ2/test_data/
    INDEX=/data/BIO3/indices/WBcel235_bwa_0.7.5a-r405/WBcel235.fa
    
    bwa mem -t16 -k 15 -T 1 $INDEX \
        ${PREFIX}/synthetic_reads.R1.fa.gz \
        ${PREFIX}/synthetic_reads.R2.fa.gz \
        > /scratch/circdetection_test_data/marcel_test.sam
    
    cat /scratch/circdetection_test_data/marcel_test.sam | \
        ../find_circ.py --test -G $INDEX --stdout=test | les
```

# Human circRNAs from circbase

```
    # in find_circ2/test_data
    # prepare a sample of human circRNAs
    wget http://www.circbase.org/download/hsa_hg19_circRNA.bed
    cat hsa_hg19_circRNA.bed | unsort.py | head -n 1000 | cut -f 1,2,3,4,5,6 | \
        exon_intersect.py -S hg19 -m hg19_sample.ucsc > hg19_sample.bed

    INDEX=/data/rajewsky/indices/hg19_bwa_0.7.12-r1039/hg19.fa
    N_FRAGS=100
    SEED=1336587
    
    # generating simulated reads
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.00 -o sim/hg19/r100_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.005 -o sim/hg19/r100_f350_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.01 -o sim/hg19/r100_f350_fpk100/m0.01 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.02 -o sim/hg19/r100_f350_fpk100/m0.02 --fpk --seed=${SEED} &

    # generating shorter simulated reads
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.00 -o sim/hg19/r50_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.005 -o sim/hg19/r50_f350_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.01 -o sim/hg19/r50_f350_fpk100/m0.01 --fpk --seed=${SEED} &

    # generating longer simulated reads
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.00 -o sim/hg19/r150_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.005 -o sim/hg19/r150_f350_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.01 -o sim/hg19/r150_f350_fpk100/m0.01 --fpk --seed=${SEED} &

    # generating very long simulated reads
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.00 -o sim/hg19/r200_f450_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.005 -o sim/hg19/r200_f450_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat hg19_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.01 -o sim/hg19/r200_f450_fpk100/m0.01 --fpk --seed=${SEED} &


    declare -a hg19_samples=(hg19/r50_f350_fpk100/m0.00 hg19/r50_f350_fpk100/m0.005 hg19/r50_f350_fpk100/m0.01 
                             hg19/r100_f350_fpk100/m0.00 hg19/r100_f350_fpk100/m0.005 hg19/r100_f350_fpk100/m0.01 hg19/r100_f350_fpk100/m0.02 
                             hg19/r150_f350_fpk100/m0.00 hg19/r150_f350_fpk100/m0.005 hg19/r150_f350_fpk100/m0.01 
                             hg19/r200_f450_fpk100/m0.00 hg19/r200_f450_fpk100/m0.005 hg19/r200_f450_fpk100/m0.01)

    # mapping
    for SAMPLE in "${hg19_samples[@]}"
    do {
        mkdir -p run/${SAMPLE}
        bwa mem -t16 -k 14 -T 1 -L 3,3 -O 6,6 -E 3,3 -p $INDEX sim/${SAMPLE}/simulated_reads.fa.gz > run/${SAMPLE}/aligned.sam
    } done;

    # running find_circ2
    for SAMPLE in "${hg19_samples[@]}"
    do {
        ../find_circ.py run/${SAMPLE}/aligned.sam -o run/${SAMPLE} --test -G $INDEX \
            --known-lin=sim/${SAMPLE}/lin_splice_sites.bed \
            --known-circ=sim/${SAMPLE}/circ_splice_sites.bed &
    } done;

    # making scatter plots
    for SAMPLE in "${hg19_samples[@]}"
    do {
        ../merge_bed.py -6 --score sim/${SAMPLE}/circ_splice_sites.bed run/${SAMPLE}/circ_splice_sites.bed | \
            histogram.py -s -q -b0 --linear-fit -x "simulated junction reads" -y "recovered junction reads" -t "backspliced read recovery" --pdf=run/${SAMPLE}/run_vs_sim.pdf
    } done;
```


# C.elegans circRNAs from circBASE
```
    wget http://www.circbase.org/download/cel_ce6_circRNA.bed
    cat cel_ce6_circRNA.bed | unsort.py | head -n 1000 | cut -f 1,2,3,4,5,6 | \
        exon_intersect.py -S ce6 -m ce6_sample.ucsc > ce6_sample.bed

    INDEX=/data/rajewsky/indices/ce6_bwa_0.7.5a/ce6.fa
    N_FRAGS=100
    SEED=1336587

    # generating simulated reads
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.00 -o sim/ce6/r100_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.005 -o sim/ce6/r100_f350_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.01 -o sim/ce6/r100_f350_fpk100/m0.01 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.02 -o sim/ce6/r100_f350_fpk100/m0.02 --fpk --seed=${SEED} &

    # generating shorter simulated reads
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.00 -o sim/ce6/r50_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.005 -o sim/ce6/r50_f350_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.01 -o sim/ce6/r50_f350_fpk100/m0.01 --fpk --seed=${SEED} &

    # generating longer simulated reads
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.00 -o sim/ce6/r150_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.005 -o sim/ce6/r150_f350_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.01 -o sim/ce6/r150_f350_fpk100/m0.01 --fpk --seed=${SEED} &

    # generating very long simulated reads
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.00 -o sim/ce6/r200_f450_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.005 -o sim/ce6/r200_f450_fpk100/m0.005 --fpk --seed=${SEED} &
    time cat ce6_sample.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.01 -o sim/ce6/r200_f450_fpk100/m0.01 --fpk --seed=${SEED} &


    declare -a ce6_samples=(ce6/r50_f350_fpk100/m0.00 ce6/r50_f350_fpk100/m0.005 ce6/r50_f350_fpk100/m0.01 
                             ce6/r100_f350_fpk100/m0.00 ce6/r100_f350_fpk100/m0.005 ce6/r100_f350_fpk100/m0.01 ce6/r100_f350_fpk100/m0.02 
                             ce6/r150_f350_fpk100/m0.00 ce6/r150_f350_fpk100/m0.005 ce6/r150_f350_fpk100/m0.01 
                             ce6/r200_f450_fpk100/m0.00 ce6/r200_f450_fpk100/m0.005 ce6/r200_f450_fpk100/m0.01)

    # mapping
    for SAMPLE in "${ce6_samples[@]}"
    do {
        mkdir -p run/${SAMPLE}
        bwa mem -t16 -k 14 -T 1 -L 3,3 -O 6,6 -E 3,3 -p $INDEX sim/${SAMPLE}/simulated_reads.fa.gz > run/${SAMPLE}/aligned.sam
    } done;

    # running find_circ2
    for SAMPLE in "${ce6_samples[@]}"
    do {
        ../find_circ.py run/${SAMPLE}/aligned.sam -o run/${SAMPLE} --test -G $INDEX \
            --known-lin=sim/${SAMPLE}/lin_splice_sites.bed \
            --known-circ=sim/${SAMPLE}/circ_splice_sites.bed &
    } done;

    # making scatter plots
    for SAMPLE in "${ce6_samples[@]}"
    do {
        ../merge_bed.py -6 --score sim/${SAMPLE}/circ_splice_sites.bed run/${SAMPLE}/circ_splice_sites.bed | \
            histogram.py -s -q -b0 --linear-fit -x "simulated junction reads" -y "recovered junction reads" -t "backspliced read recovery" --pdf=run/${SAMPLE}/run_vs_sim.pdf
    } done;

```



# C.elegans 'golden set'

Sensitivity seems to decrease in C.elegans when read length is 200. This is quite strange and we suspect it could be due to false positives in circBASE. For internal use, we generate a 'golden set' of circ-RNAs that don't get any warnings by find_circ2.

```
    INDEX=/data/rajewsky/indices/ce6_bwa_0.7.5a/ce6.fa
    N_FRAGS=100
    SEED=1336587

    cat run/ce6/r200_f450_fpk100/m0.00/circ_splice_sites.bed | grep -v WARN_EXT_2MM | grep -v WARN_OUTSIDE | grep -v WARN_MULTI > ce6_golden_set.bed
    cut -f 1,2,3,4,5,6 ce6_golden_set.bed | grep -v '#' | exon_intersect.py -S ce6 -m ce6_golden.ucsc > ce6_golden.bed
    
    # generating simulated reads
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.00 -o sim_gold/ce6/r100_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.01 -o sim_gold/ce6/r100_f350_fpk100/m0.01 --fpk --seed=${SEED} &
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.05 -o sim_gold/ce6/r100_f350_fpk100/m0.05 --fpk --seed=${SEED} &
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.10 -o sim_gold/ce6/r100_f350_fpk100/m0.10 --fpk --seed=${SEED} &

    # generating shorter simulated reads
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.00 -o sim_gold/ce6/r50_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.01 -o sim_gold/ce6/r50_f350_fpk100/m0.01 --fpk --seed=${SEED} &

    # generating longer simulated reads
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.00 -o sim_gold/ce6/r150_f350_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.01 -o sim_gold/ce6/r150_f350_fpk100/m0.01 --fpk --seed=${SEED} &

    # generating very long simulated reads
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.00 -o sim_gold/ce6/r200_f450_fpk100/m0.00 --fpk --seed=${SEED} &
    time cat ce6_golden.ucsc | grep ANNOTATED | grep INTERNAL | \
        ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.01 -o sim_gold/ce6/r200_f450_fpk100/m0.01 --fpk --seed=${SEED} &


    declare -a ce6_samples=(ce6/r50_f350_fpk100/m0.01 ce6/r100_f350_fpk100/m0.00 ce6/r100_f350_fpk100/m0.01 ce6/r100_f350_fpk100/m0.05 ce6/r100_f350_fpk100/m0.10 ce6/r150_f350_fpk100/m0.01 ce6/r200_f450_fpk100/m0.01 ce6/r50_f350_fpk100/m0.00 ce6/r150_f350_fpk100/m0.00 ce6/r200_f450_fpk100/m0.00)

    # mapping
    for SAMPLE in "${ce6_samples[@]}"
    do {
        mkdir -p run_gold/${SAMPLE}
        bwa mem -t16 -k 14 -T 1 -L 3,3 -O 6,6 -E 3,3 -p $INDEX sim_gold/${SAMPLE}/simulated_reads.fa.gz > run_gold/${SAMPLE}/aligned.sam
    } done;

    # running find_circ2
    for SAMPLE in "${ce6_samples[@]}"
    do {
        ../find_circ.py run_gold/${SAMPLE}/aligned.sam -o run_gold/${SAMPLE} --test -G $INDEX \
            --known-lin=sim/${SAMPLE}/lin_splice_sites.bed \
            --known-circ=sim/${SAMPLE}/circ_splice_sites.bed &
    } done;

    # making scatter plots
    for SAMPLE in "${ce6_samples[@]}"
    do {
        ../merge_bed.py -6 --score sim_gold/${SAMPLE}/circ_splice_sites.bed run_gold/${SAMPLE}/circ_splice_sites.bed | \
            histogram.py -s -q -b0 --linear-fit -x "simulated junction reads" -y "recovered junction reads" -t "backspliced read recovery" --pdf=run_gold/${SAMPLE}/run_vs_sim.pdf
    } done;

```


# let's find out which circRNAs have higher/lower sensitivity of circRNAs of different lengths
```
    declare -a ce6_samples=(ce6/r50_f350_fpk100/m0.00 ce6/r100_f350_fpk100/m0.00 ce6/r150_f350_fpk100/m0.00 ce6/r200_f450_fpk100/m0.00)
    for SAMPLE in "${ce6_samples[@]}"
    do {
        ../merge_bed.py -6 --score sim_gold/${SAMPLE}/circ_splice_sites.bed run_gold/${SAMPLE}/circ_splice_sites.bed | \
            lfc.py 2 3 1 | cut -f 1,5 > run_gold/${SAMPLE}/circ_sens.tsv
    } done;
    
    merge.py run_gold/ce6/r150_f350_fpk100/m0.00/circ_sens.tsv run_gold/ce6/r200_f450_fpk100/m0.00/circ_sens.tsv -u | \
        histogram.py -s -x "recovery r150" -y "recovery r200" -S -q --pdf funky_recovery.pdf
    
    merge.py run_gold/ce6/r150_f350_fpk100/m0.00/circ_sens.tsv run_gold/ce6/r200_f450_fpk100/m0.00/circ_sens.tsv -u | \
        scorethresh.py 2 0 | scorethresh.py -3 -.5 > two_funky_circs
```