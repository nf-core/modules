process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2cb827988dca7067ff8096c37cb20bc841c878013da52ad47a50865d54efe83/data' :
        'community.wave.seqera.io/library/megahit_pigz:87a590163e594224' }"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    tuple val(meta), path("intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    tuple val(meta), path("intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    tuple val(meta), path("intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    tuple val(meta), path('*.log')                                      , emit: log
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end || !reads2 ? "-r ${reads1.join(',')}" : "-1 ${reads1.join(',')} -2 ${reads2.join(',')}"
    """
    megahit \\
        ${args} \\
        -t ${task.cpus} \\
        ${reads_command} \\
        --out-prefix ${prefix}

    pigz \\
        --no-name \\
        -p ${task.cpus} \\
        ${args2} \\
        megahit_out/*.fa \\
        megahit_out/intermediate_contigs/*.fa

    mv megahit_out/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end || !reads2 ? "-r ${reads1}" : "-1 ${reads1.join(',')} -2 ${reads2.join(',')}"
    """
    mkdir -p intermediate_contigs
    echo "" | gzip > ${prefix}.contigs.fa.gz
    echo "" | gzip > intermediate_contigs/k21.contigs.fa.gz
    echo "" | gzip > intermediate_contigs/k21.addi.fa.gz
    echo "" | gzip > intermediate_contigs/k21.local.fa.gz
    echo "" | gzip > intermediate_contigs/k21.final.contigs.fa.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
