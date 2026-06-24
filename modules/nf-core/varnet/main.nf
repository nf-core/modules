process VARNET {
    tag "$meta.id"
    label 'process_high_memory'

    container "docker.io/kiranchari/varnet:latest"

    input:
    tuple val(meta), path(input_tumor), path(index_tumor), path(input_normal), path(index_normal)
    tuple val(meta2), path(intervals)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("${prefix}/${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val("varnet"), env("VARNET_VERSION"), emit: versions_varnet, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--region_bed \$WORKDIR/${intervals}" : ""
    def normal  = input_normal ? "--normal_bam \$WORKDIR/${input_normal}" : ""
    """

    WORKDIR=\$(pwd)

    cd /VarNet
    export VARNET_VERSION=\$(python -c "from snvs.constants import __VERSION__; print(__VERSION__)")
    TF_CPP_MIN_LOG_LEVEL=3 python /VarNet/filter.py \\
        --sample_name ${prefix} \\
        ${normal} \\
        --tumor_bam \$WORKDIR/${input_tumor} \\
        --reference \$WORKDIR/${fasta} \\
        --output_dir \$WORKDIR \\
        --processes ${task.cpus} \\
        ${regions} \\
        ${args}

    TF_CPP_MIN_LOG_LEVEL=3 python /VarNet/predict.py \\
        --sample_name ${prefix} \\
        ${normal} \\
        --tumor_bam \$WORKDIR/${input_tumor} \\
        --reference \$WORKDIR/${fasta} \\
        --output_dir \$WORKDIR \\
        --processes ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export VARNET_VERSION=\$(python -c "from snvs.constants import __VERSION__; print(__VERSION__)")
    mkdir -p ${prefix}
    echo "" | gzip > ${prefix}/${prefix}.vcf.gz
    """
}

