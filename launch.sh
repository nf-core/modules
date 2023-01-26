#!/bin/bash
conda activate env_nf
nohup nextflow run ./tests/modules/nf-core/glimpse/ligate -entry test_glimpse_ligate -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/glimpse/ligate/nextflow.config --outdir ./data -work-dir ./work --resume &> nohupBQSR.out