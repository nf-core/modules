nextflow run ~/repos/nf-core/modules/tests/modules/nf-core/somatosim/somatosim_workflow.nf -with-docker somatosim:latest -dump-channels -resume -dump-hashes

somatosim -i ${prefix}.bam \\
        -b ${bed} \\
        -o output_dir \\
        --vaf-low 0.01 \\
        --vaf-high 0.05 \\
        --number-snv 100 \\
        --random-seed 0 \\
        --verbose