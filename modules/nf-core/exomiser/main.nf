process EXOMISER {
    tag "$meta.id"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::exomiser-cli=13.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/exomiser-cli:13.2.1--hdfd78af_0':
        'biocontainers/exomiser-cli:13.2.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), assembly, path(analysis_yml), path(job_yml) , path(output_yml), path(ped), path(sample_phenopacket_yml)

    output:
    tuple val(meta), path("*.vcf.gz")       , optional:true, emit: vcf
    tuple val(meta), path("*.html")         , optional:true, emit: html
    tuple val(meta), path("*.json")         , optional:true, emit: json
    tuple val(meta), path("*genes.tsv")     , optional:true, emit: genes_tsv
    tuple val(meta), path("*variants.tsv")  , optional:true, emit: variants_tsv

    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def EXOMISER_VERSION = 13.2.1
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    exomiser-cli \\
        --analysis ${analysis_yml} \\
        --assembly ${assembly} \\
        --job ${job_yml} \\
        --output ${output_yml} \\
        --output-directory . \\
        --output-filename ${prefix} \\
        --ped ${ped} \\
        --sample ${sample_phenopacket_yml} \\
        --vcf ${vcf} \\
        $args

    ## $args can contain
    ## --output-format
    ## --output-prefix (deprecated)
    ## --preset (mandatory)


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : ${EXOMISER_VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf
    touch ${prefix}.html
    touch ${prefix}.json
    touch ${prefix}.genes.tsv
    touch ${prefix}.variants.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : ${EXOMISER_VERSION}
    END_VERSIONS
    """
}
