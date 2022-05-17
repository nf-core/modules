#!/usr/bin/env bash
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
# Copyright (c) 2022 Genome Research Ltd.
#
# =============================================================================
# Setup Instructions
# =============================================================================
#
#    Add any other necessary LSF arguments such as queue (-q) or account (-P).
#    If your system requires a walltime (-W), 24 hours (24:00) is sufficient.
#    We recommend you do not remove any arguments below or Martian may not run
#    properly.
#
# =============================================================================
# Template
# =============================================================================
#
#BSUB -J __MRO_JOB_NAME__
#BSUB -n __MRO_THREADS__
#BSUB -o __MRO_STDOUT__
#BSUB -e __MRO_STDERR__
#BSUB -R "select[type==X86_64 && mem>__MRO_MEM_MB__] rusage[mem=__MRO_MEM_MB__]"
#BSUB -R span[hosts=1]
#BSUB -q long
#BSUB -M __MRO_MEM_MB__

singularity exec ${SINGULARITY_CONTAINER} env PATH="$PATH" \
__MRO_CMD__
