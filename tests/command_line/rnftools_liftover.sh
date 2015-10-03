#! /usr/bin/env bash

set -eux
set -o pipefail

cd "$(dirname "$0")"

rnftools liftover --genome-id 1 -c liftover.chain liftover_orig.bam - > _liftover_new1.sam
rnftools liftover --genome-id 1 -c liftover.chain --output-format bam liftover_orig.bam - > _liftover_new1.bam

rnftools liftover --genome-id 1 liftover_orig.bam _liftover_orig.sam

rnftools liftover --genome-id 1 -c liftover.chain liftover_orig.bam _liftover_new2.sam
rnftools liftover --genome-id 1 -c liftover.chain liftover_orig.bam _liftover_new2.bam
rnftools liftover --genome-id 1 -c liftover.chain _liftover_orig.sam _liftover_new3.sam
rnftools liftover --genome-id 1 -c liftover.chain _liftover_orig.sam _liftover_new3.bam

rnftools liftover --genome-id 1 --invert -c liftover.chain _liftover_new1.bam _liftover_pseudoorig.sam
