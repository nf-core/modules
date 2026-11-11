process HGTSITE_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hgt-site:0.4.0--pyhdfd78af_0' :
        'quay.io/biocontainers/hgt-site:0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(alignment), path(index)
    tuple val(meta2), path(classifications)
    val host_taxid
    val save_bed
    val save_vcf
    val save_multiqc

    output:
    tuple val(meta), path("*.tsv")                            , emit: tsv
    tuple val(meta), path("*.bed")            , optional:true , emit: bed
    tuple val(meta), path("*.vcf.gz")         , optional:true , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi")     , optional:true , emit: tbi
    tuple val(meta), path("*_mqc.json")       , optional:true , emit: multiqc
    tuple val("${task.process}"), val('hgtsite'), eval('hgt-site --version | sed "s/^hgt-site //"'), topic: versions, emit: versions_hgtsite

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_option = save_bed ? "--bed ${prefix}.bed" : ""
    def vcf_option = save_vcf ? "--vcf ${prefix}.vcf.gz" : ""
    def multiqc_option = save_multiqc ? "--multiqc ${prefix}_mqc.json" : ""
    """
    hgt-site call \\
        --bam ${alignment} \\
        --classifications ${classifications} \\
        --host-taxid ${host_taxid} \\
        --threads ${task.cpus} \\
        --sample-name ${prefix} \\
        ${bed_option} \\
        ${vcf_option} \\
        ${multiqc_option} \\
        -o ${prefix}.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed_stub = save_bed ? "touch ${prefix}.bed" : ""
    def vcf_stub = save_vcf ? "echo '' | gzip > ${prefix}.vcf.gz && touch ${prefix}.vcf.gz.tbi" : ""
    def multiqc_stub = save_multiqc ? "touch ${prefix}_mqc.json" : ""
    """
    touch ${prefix}.tsv
    ${bed_stub}
    ${vcf_stub}
    ${multiqc_stub}
    """
}
