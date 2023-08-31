process EXOMISER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::exomiser-rest-prioritiser=13.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/exomiser-rest-prioritiser:13.2.0--hdfd78af_0':
        'biocontainers/exomiser-rest-prioritiser:13.2.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(config)
    tuple val(meta) , path(phenopacket)
    tuple val(meta) , path(analysis_parameters)
    tuple val(meta) , path(output_parameters)

    output:
    tuple val(meta), path("*.vcf")          , optional:true, emit: vcf
    tuple val(meta), path("*.html")         , optional:true, emit: html
    tuple val(meta), path("*.json")         , optional:true, emit: json
    tuple val(meta), path("*.vcf.gz.tbi")   , optional:true, emit: vcf.gz.tbi
    tuple val(meta), path("*genes.tsv")     , optional:true, emit: genes.tsv
    tuple val(meta), path("*variants.tsv")  , optional:true, emit: variants.tsv

    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Normal command: java -jar exomiser-cli-13.2.0.jar --sample examples/pfeiffer-phenopacket.yml --analysis examples/exome-analysis.yml --output examples/output-options.yml
    # @Matthias: Would sample and analysis flag be somehow handled by $args?

    exomiser/exomiser-cli:13.2.0 \\
        $args \\
        --sample ${phenopacket} \\
        --analysis ${analysis_parameters} \\
        --output ${output_parameters}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(exomiser --version 2>&1) | sed 's/^.*exomiser //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(exomiser --version 2>&1) | sed 's/^.*exomiser //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
