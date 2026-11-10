process HGTSITE_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hgt-site:0.4.0--pyhdfd78af_0' :
        'biocontainers/hgt-site:0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(alignment), path(index)
    tuple val(meta2), path(classifications)
    val host_taxid

    output:
    tuple val(meta), path("*.hgt_site.tsv")                                     , emit: tsv
    tuple val(meta), path("*.hgt_site.bed")                      , optional:true, emit: bed
    tuple val(meta), path("*.hgt_site.vcf{,.gz}")                , optional:true, emit: vcf
    tuple val(meta), path("*.hgt_site.vcf.gz.tbi")               , optional:true, emit: tbi
    tuple val(meta), path("*_mqc.json")                          , optional:true, emit: multiqc
    tuple val("${task.process}"), val('hgtsite'), eval('hgt-site --version | sed "s/^hgt-site //"'), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    hgt-site call \\
        --bam ${alignment} \\
        --classifications ${classifications} \\
        --host-taxid ${host_taxid} \\
        --threads ${task.cpus} \\
        --sample-name ${prefix} \\
        -o ${prefix}.hgt_site.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hgt_site.tsv
    """
}
