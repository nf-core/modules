process TOPAS_GENCONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/topas:1.0.1--hdfd78af_1':
        'biocontainers/topas:1.0.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(vcf_indels)
    tuple val(meta3), path(reference)
    tuple val(meta4), path(fai)
    val(vcf_output)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.vcf.gz")  , emit: vcf     , optional: true
    tuple val(meta), path("*.ccf")     , emit: ccf
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def optionalvcfindels = vcf_indels ? "-indels ${vcf_indels}" : ''
    def optionalfai = fai ? "-fai ${fai}" : ''
    def vcfoutput = vcf_output ? "-vcf_out ${prefix}.vcf" : ""
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """

    topas \\
        GenConS \\
        $args \\
        -o ${prefix}.fasta \\
        -snps $vcf \\
        $optionalvcfindels \\
        $optionalfai \\
        $vcfoutput \\
        -ref $reference

    gzip -n ${prefix}.fasta

    if [[ -f ${prefix}.vcf ]];then
        gzip -n ${prefix}.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        topas: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcfoutput = vcf_output ? "echo | gzip > ${prefix}.vcf.gz" : ""
    def VERSION = '1.0.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo | gzip > ${prefix}.fasta.gz
    touch ${prefix}.fastq.ccf
    touch ${prefix}.fastq.log
    $vcfoutput

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        topas: $VERSION
    END_VERSIONS
    """

}
