#!/bin/bash
set -e -x

kalign -q -f clu -s 100 -e 0.85 -t 0.45 -m 0 -in ${dir}/seq.fa \
    | sed 's/Kalign/CLUSTAL/g' > ${dir}/msa.clu

python ${python_scripts}/msa_to_variants.py \
    ${dir}/msa.clu \
    ${chrom} \
    ${start} \
    ${end} \
    ${dir}/seq.qual \
    50 > ${dir}/svs.bed

