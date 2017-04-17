#! /usr/bin/env bash

set -euf -o pipefail

cd "$(dirname "$0")"

# single end
rnftools merge \
	-i rnf_fq_se_1.fq rnf_fq_se_2.fq \
	-m single-end \
	-o _mix_se

# paired-end
rnftools merge \
	-i rnf_fq_pe_1.fq rnf_fq_pe_2.fq \
	-m paired-end-bwa \
	-o _mix_pe

rnftools merge \
	-i rnf_fq_pe_1.fq rnf_fq_pe_2.fq \
	-m paired-end-bfast \
	-o _mix_pe