# Measuring find_circ2 performance

This worksheet describes the simulations and analyses run to make the figures for *Jens et al. 2016* .

## Testing sensitivity and precision

First we simulate reads from known circRNAs in human and worm to benchmark the sensitivity and precision of spliced read recovery.

### Simulating human circRNAs from circbase

We take 1000 circRNAs randomly from circBase. Analysis is further restricted to circRNAs of at least 200nt to have a comparable set from which different read-lengths can be drawn. The simulation does currently not allow drawing reads that are longer than the circRNA.

```
    # in find_circ2/test_data
    # prepare a sample of human circRNAs
    wget http://www.circbase.org/download/hsa_hg19_circRNA.bed
    cat hsa_hg19_circRNA.bed | unsort.py | head -n 1000 | cut -f 1,2,3,4,5,6 | \
        exon_intersect.py -S hg19 -m hg19_sample.ucsc > hg19_sample.bed

    # get 1,000 circRNA exon counts for the generation of random circ structures.
    cat hg19_sample.ucsc | grep -v '#' | head -n 1000 | cut -f 9 | sort -grk 1 > hg19.circ.exon_counts
    histogram.py -b 0 hg19.circ.exon_counts -p hg19.circ.exon_counts.pdf

    # obtained GENCODEV19 from UCSC table browser.
    zcat gencodev19.gtf.gz | grep exon > hg19.gencode_v19.exons.gtf
    
    for SEED in 1337 4711 815 110112 11358
    do {
        cat hg19.circ.exon_counts | \
            python ../random_circs.py hg19.gencode_v19.exons.gtf --seed ${SEED} \
                > hg19.synth_circs.${SEED}.ucsc \
                2> hg19.synth_circs.${SEED}.log &
    } done;

    
    INDEX=/data/rajewsky/indices/hg19_bwa_0.7.12-r1039/hg19.fa
    N_FRAGS=100
    SEED=1336587
    
    for SEED in 1337 4711 815 110112 11358
    do {  
        # generating simulated reads
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.00 -o sim/hg19/${SEED}/r100_f350_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.005 -o sim/hg19/${SEED}/r100_f350_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.01 -o sim/hg19/${SEED}/r100_f350_fpk100/m0.01 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.02 -o sim/hg19/${SEED}/r100_f350_fpk100/m0.02 --fpk --seed=${SEED} &

        # generating shorter simulated reads
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.00 -o sim/hg19/${SEED}/r50_f350_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.005 -o sim/hg19/${SEED}/r50_f350_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.01 -o sim/hg19/${SEED}/r50_f350_fpk100/m0.01 --fpk --seed=${SEED} &

        # generating longer simulated reads
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.00 -o sim/hg19/${SEED}/r150_f350_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.005 -o sim/hg19/${SEED}/r150_f350_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.01 -o sim/hg19/${SEED}/r150_f350_fpk100/m0.01 --fpk --seed=${SEED} &

        # generating very long simulated reads
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.00 -o sim/hg19/${SEED}/r200_f450_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.005 -o sim/hg19/${SEED}/r200_f450_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat hg19.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.01 -o sim/hg19/${SEED}/r200_f450_fpk100/m0.01 --fpk --seed=${SEED}
    } done;

    declare -a hg19_samples=(
        hg19/1337/r50_f350_fpk100/m0.00 hg19/1337/r50_f350_fpk100/m0.005 hg19/1337/r50_f350_fpk100/m0.01 
        hg19/1337/r100_f350_fpk100/m0.00 hg19/1337/r100_f350_fpk100/m0.005 hg19/1337/r100_f350_fpk100/m0.01 hg19/1337/r100_f350_fpk100/m0.02 
        hg19/1337/r150_f350_fpk100/m0.00 hg19/1337/r150_f350_fpk100/m0.005 hg19/1337/r150_f350_fpk100/m0.01 
        hg19/1337/r200_f450_fpk100/m0.00 hg19/1337/r200_f450_fpk100/m0.005 hg19/1337/r200_f450_fpk100/m0.01
        hg19/4711/r50_f350_fpk100/m0.00 hg19/4711/r50_f350_fpk100/m0.005 hg19/4711/r50_f350_fpk100/m0.01 
        hg19/4711/r100_f350_fpk100/m0.00 hg19/4711/r100_f350_fpk100/m0.005 hg19/4711/r100_f350_fpk100/m0.01 hg19/4711/r100_f350_fpk100/m0.02 
        hg19/4711/r150_f350_fpk100/m0.00 hg19/4711/r150_f350_fpk100/m0.005 hg19/4711/r150_f350_fpk100/m0.01 
        hg19/4711/r200_f450_fpk100/m0.00 hg19/4711/r200_f450_fpk100/m0.005 hg19/4711/r200_f450_fpk100/m0.01
        hg19/815/r50_f350_fpk100/m0.00 hg19/815/r50_f350_fpk100/m0.005 hg19/815/r50_f350_fpk100/m0.01 
        hg19/815/r100_f350_fpk100/m0.00 hg19/815/r100_f350_fpk100/m0.005 hg19/815/r100_f350_fpk100/m0.01 hg19/815/r100_f350_fpk100/m0.02 
        hg19/815/r150_f350_fpk100/m0.00 hg19/815/r150_f350_fpk100/m0.005 hg19/815/r150_f350_fpk100/m0.01 
        hg19/815/r200_f450_fpk100/m0.00 hg19/815/r200_f450_fpk100/m0.005 hg19/815/r200_f450_fpk100/m0.01
        hg19/110112/r50_f350_fpk100/m0.00 hg19/110112/r50_f350_fpk100/m0.005 hg19/110112/r50_f350_fpk100/m0.01 
        hg19/110112/r100_f350_fpk100/m0.00 hg19/110112/r100_f350_fpk100/m0.005 hg19/110112/r100_f350_fpk100/m0.01 hg19/110112/r100_f350_fpk100/m0.02 
        hg19/110112/r150_f350_fpk100/m0.00 hg19/110112/r150_f350_fpk100/m0.005 hg19/110112/r150_f350_fpk100/m0.01 
        hg19/110112/r200_f450_fpk100/m0.00 hg19/110112/r200_f450_fpk100/m0.005 hg19/110112/r200_f450_fpk100/m0.01
        hg19/11358/r50_f350_fpk100/m0.00 hg19/11358/r50_f350_fpk100/m0.005 hg19/11358/r50_f350_fpk100/m0.01 
        hg19/11358/r100_f350_fpk100/m0.00 hg19/11358/r100_f350_fpk100/m0.005 hg19/11358/r100_f350_fpk100/m0.01 hg19/11358/r100_f350_fpk100/m0.02 
        hg19/11358/r150_f350_fpk100/m0.00 hg19/11358/r150_f350_fpk100/m0.005 hg19/11358/r150_f350_fpk100/m0.01 
        hg19/11358/r200_f450_fpk100/m0.00 hg19/11358/r200_f450_fpk100/m0.005 hg19/11358/r200_f450_fpk100/m0.01
    )

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
            --known-circ=sim/${SAMPLE}/circ_splice_sites.bed 2> run/${SAMPLE}/jk2 &
    } done;

    # making scatter plots
    for SAMPLE in "${hg19_samples[@]}"
    do {
        ../merge_bed.py -6 --score sim/${SAMPLE}/circ_splice_sites.bed run/${SAMPLE}/circ_splice_sites.bed | \
            histogram.py -s -q -b0 --linear-fit -x "simulated junction reads" -y "recovered junction reads" -t "backspliced read recovery" --pdf=run/${SAMPLE}/run_vs_sim.pdf
    } done;
```


### C.elegans circRNAs from circBASE

We repeat this exercise for *C.elegans*

```
    wget http://www.circbase.org/download/cel_ce6_circRNA.bed
    cat cel_ce6_circRNA.bed | unsort.py | head -n 1000 | cut -f 1,2,3,4,5,6 | \
        exon_intersect.py -S ce6 -m ce6_sample.ucsc > ce6_sample.bed

    # get 1,000 circRNA exon counts for the generation of random circ structures.
    cat ce6_sample.ucsc ce6_sample.ucsc | unsort.py | grep -v '#' | head -n 1000 | cut -f 9 | sort -grk 1 > ce6.circ.exon_counts
    histogram.py -b 0 ce6.circ.exon_counts -p ce6.circ.exon_counts.pdf

    
    for SEED in 1337 4711 815 110112 11358
    do {
        cat ce6.circ.exon_counts | \
            python ../random_circs.py ce6.ensGene.exons.gtf --seed ${SEED} \
                > ce6.synth_circs.${SEED}.ucsc \
                2> ce6.synth_circs.${SEED}.log &
    } done;
    
    INDEX=/data/rajewsky/indices/ce6_bwa_0.7.5a/ce6.fa
    N_FRAGS=100

    for SEED in 1337 4711 815 110112 11358
    do {  
        # generating simulated reads
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.00 -o sim/ce6/${SEED}/r100_f350_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.005 -o sim/ce6/${SEED}/r100_f350_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.01 -o sim/ce6/${SEED}/r100_f350_fpk100/m0.01 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=100 --frag-len=350 --mutate=0.02 -o sim/ce6/${SEED}/r100_f350_fpk100/m0.02 --fpk --seed=${SEED} &

        # generating shorter simulated reads
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.00 -o sim/ce6/${SEED}/r50_f350_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.005 -o sim/ce6/${SEED}/r50_f350_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=50 --frag-len=350 --mutate=0.01 -o sim/ce6/${SEED}/r50_f350_fpk100/m0.01 --fpk --seed=${SEED} &

        # generating longer simulated reads
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.00 -o sim/ce6/${SEED}/r150_f350_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.005 -o sim/ce6/${SEED}/r150_f350_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=150 --frag-len=350 --mutate=0.01 -o sim/ce6/${SEED}/r150_f350_fpk100/m0.01 --fpk --seed=${SEED} &

        # generating very long simulated reads
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.00 -o sim/ce6/${SEED}/r200_f450_fpk100/m0.00 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.005 -o sim/ce6/${SEED}/r200_f450_fpk100/m0.005 --fpk --seed=${SEED} &
        time cat ce6.synth_circs.${SEED}.ucsc | \
            ../simulate_reads.py -G $INDEX --n-frags=${N_FRAGS} --read-len=200 --frag-len=450 --mutate=0.01 -o sim/ce6/${SEED}/r200_f450_fpk100/m0.01 --fpk --seed=${SEED}
    } done;

    declare -a ce6_samples=(
        ce6/1337/r50_f350_fpk100/m0.00 ce6/1337/r50_f350_fpk100/m0.005 ce6/1337/r50_f350_fpk100/m0.01 
        ce6/1337/r100_f350_fpk100/m0.00 ce6/1337/r100_f350_fpk100/m0.005 ce6/1337/r100_f350_fpk100/m0.01 ce6/1337/r100_f350_fpk100/m0.02 
        ce6/1337/r150_f350_fpk100/m0.00 ce6/1337/r150_f350_fpk100/m0.005 ce6/1337/r150_f350_fpk100/m0.01 
        ce6/1337/r200_f450_fpk100/m0.00 ce6/1337/r200_f450_fpk100/m0.005 ce6/1337/r200_f450_fpk100/m0.01
        ce6/4711/r50_f350_fpk100/m0.00 ce6/4711/r50_f350_fpk100/m0.005 ce6/4711/r50_f350_fpk100/m0.01 
        ce6/4711/r100_f350_fpk100/m0.00 ce6/4711/r100_f350_fpk100/m0.005 ce6/4711/r100_f350_fpk100/m0.01 ce6/4711/r100_f350_fpk100/m0.02 
        ce6/4711/r150_f350_fpk100/m0.00 ce6/4711/r150_f350_fpk100/m0.005 ce6/4711/r150_f350_fpk100/m0.01 
        ce6/4711/r200_f450_fpk100/m0.00 ce6/4711/r200_f450_fpk100/m0.005 ce6/4711/r200_f450_fpk100/m0.01
        ce6/815/r50_f350_fpk100/m0.00 ce6/815/r50_f350_fpk100/m0.005 ce6/815/r50_f350_fpk100/m0.01 
        ce6/815/r100_f350_fpk100/m0.00 ce6/815/r100_f350_fpk100/m0.005 ce6/815/r100_f350_fpk100/m0.01 ce6/815/r100_f350_fpk100/m0.02 
        ce6/815/r150_f350_fpk100/m0.00 ce6/815/r150_f350_fpk100/m0.005 ce6/815/r150_f350_fpk100/m0.01 
        ce6/815/r200_f450_fpk100/m0.00 ce6/815/r200_f450_fpk100/m0.005 ce6/815/r200_f450_fpk100/m0.01
        ce6/110112/r50_f350_fpk100/m0.00 ce6/110112/r50_f350_fpk100/m0.005 ce6/110112/r50_f350_fpk100/m0.01 
        ce6/110112/r100_f350_fpk100/m0.00 ce6/110112/r100_f350_fpk100/m0.005 ce6/110112/r100_f350_fpk100/m0.01 ce6/110112/r100_f350_fpk100/m0.02 
        ce6/110112/r150_f350_fpk100/m0.00 ce6/110112/r150_f350_fpk100/m0.005 ce6/110112/r150_f350_fpk100/m0.01 
        ce6/110112/r200_f450_fpk100/m0.00 ce6/110112/r200_f450_fpk100/m0.005 ce6/110112/r200_f450_fpk100/m0.01
        ce6/11358/r50_f350_fpk100/m0.00 ce6/11358/r50_f350_fpk100/m0.005 ce6/11358/r50_f350_fpk100/m0.01 
        ce6/11358/r100_f350_fpk100/m0.00 ce6/11358/r100_f350_fpk100/m0.005 ce6/11358/r100_f350_fpk100/m0.01 ce6/11358/r100_f350_fpk100/m0.02 
        ce6/11358/r150_f350_fpk100/m0.00 ce6/11358/r150_f350_fpk100/m0.005 ce6/11358/r150_f350_fpk100/m0.01 
        ce6/11358/r200_f450_fpk100/m0.00 ce6/11358/r200_f450_fpk100/m0.005 ce6/11358/r200_f450_fpk100/m0.01
    )

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

## Reconstructing circRNA structure

The multievent output allows to attempt a reconstrucion of the multi-exonic structure of each circRNA. We use the simulated circRNAs as a benchmark, starting from all known exons, or no known exons at all.

```
    # obtained the ensGene transcript models from UCSC table browser in GTF format
    grep exon ce6.ensGene.gtf > ce6.ensGene.exons.gtf
    
    # de novo reconstruction
    for SAMPLE in "${ce6_samples[@]}"
    do {
            ../reconstruct_isoforms.py run/${SAMPLE}/multi_events.tsv -o run/${SAMPLE} &
    } done;

    for SAMPLE in "${hg19_samples[@]}"
    do {
            ../reconstruct_isoforms.py run/${SAMPLE}/multi_events.tsv -o run/${SAMPLE} &
    } done;
    
    # with prior knowledge of annotated exons
    ../reconstruct_isoforms.py run/ce6/r200_f450_fpk100/m0.00/multi_events.tsv --known-exons=ce6.ensGene.exons.gtf
    
```
