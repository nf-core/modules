// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CADD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/clinicalgenomics/cadd-with-scripts:1.7.3'

    containerOptions {
        if (prescored_dir) {
            ['singularity', 'apptainer'].contains(workflow.containerEngine)
                ? "-B ${annotation_dir}:/cadd-scripts/data/annotations -B ${prescored_dir}:/cadd-scripts/data/prescored"
                : "-v ${annotation_dir}:/cadd-scripts/data/annotations -v ${prescored_dir}:/cadd-scripts/data/prescored"
        }
        else {
            ['singularity', 'apptainer'].contains(workflow.containerEngine)
                ? "-B ${annotation_dir}:/cadd-scripts/data/annotations"
                : "-v ${annotation_dir}:/cadd-scripts/data/annotations"
        }
    }

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), val(annotation_dir)
    tuple val(meta3), val(prescored_dir)

    output:
    tuple val(meta), path("${prefix}.tsv.gz"), emit: tsv
    tuple val("${task.process}"), val("cadd"), val("1.7.3"), emit: versions_cadd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME=\$PWD/snakemake_cache
    export MPLCONFIGDIR=.
    mkdir -p \$XDG_CACHE_HOME

    CADD.sh \\
        -m \\
        -o ${prefix}.tsv.gz \\
        ${args} \\
        ${vcf}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.tsv.gz
    """
}
