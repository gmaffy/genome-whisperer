#!/bin/bash

for s in ON1708 ON2086; do
  L=~/gw-work/batch/data/onion/2019/$s/reference_genomes/GCA_030765085.1_ASM3076508v1/bams
  N=/mnt/v/DATA/onion/2019/$s/reference_genomes/GCA_030765085.1_ASM3076508v1/bams
  rsync -rt -h --partial --info=progress2 \
    --include='*.cram' --include='*.crai' --include='*.txt' --include='*.pdf' --exclude='*' \
    "$L"/ "$N"/
done
