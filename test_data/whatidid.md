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

# Using my own dm6 reads

```
    # in find_circ2/test_data
    INDEX=/data/rajewsky/indices/dm6_bwa_0.7.12-r1039/dm6.fa
    
    time cat ~/circpeptides/fly/dm6.ribocircs_mar2016.ucsc | grep ANNOTATED | \
        ../simulate_reads.py -G $INDEX -o sim_dm6 --fpk &

    time cat ~/circpeptides/fly/dm6.ribocircs_mar2016.ucsc | grep ANNOTATED | \
        ../simulate_reads.py -G $INDEX -o sim_dm6_0.01 --mutate=0.01 --fpk &

    time cat ~/circpeptides/fly/dm6.ribocircs_mar2016.ucsc | grep ANNOTATED | \
        ../simulate_reads.py -G $INDEX -o sim_dm6_0.05 --mutate=0.05 --fpk &

    time cat ~/circpeptides/fly/dm6.ribocircs_mar2016.ucsc | grep ANNOTATED | \
        ../simulate_reads.py -G $INDEX -o sim_dm6_0.1 --mutate=0.1 --fpk &


    for SAMPLE in sim_dm6 sim_dm6_0.01 sim_dm6_0.05 sim_dm6_0.1
    do {
        bwa mem -t16 -k 15 -T 1 -p $INDEX ${SAMPLE}/simulated_reads.fa.gz > ${SAMPLE}.sam
    
        ../find_circ.py ${SAMPLE}.sam -o ${SAMPLE}_run --test -G $INDEX \
            --known-lin=${SAMPLE}/lin_splice_sites.bed \
            --known-circ=${SAMPLE}/circ_splice_sites.bed
    } done;
    
    for ERR in "" _0.01 _0.05 _0.1
    do {
        ../merge_bed.py -6 -V sim_dm6${ERR}/circ_splice_sites.bed sim_dm6${ERR}_run/circ_splice_sites.bed | grep '(in0,in1)' | cut -f 7,14 | histogram.py -s -q -b0 --ofs-fit -x "simulated junction reads" -y "recovered junction reads" -t "circRNA recovery error=${ERR}" --pdf=scatter_sim${ERR}.pdf
    } done;
    
```