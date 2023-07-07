process ARCASHLA_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::arcas-hla=0.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arcas-hla:0.5.0--hdfd78af_1':
        'biocontainers/arcas-hla:0.5.0--hdfd78af_1' }"
    
    containerOptions = "--user root"

    input:
    tuple val(meta), path(fq)

    output:
    //tuple val(meta), path("${meta.id}/*")                     , emit: results
    // tuple val(meta), path("${meta.id}/NA12878.genes.json")        , emit: genes
    // tuple val(meta), path("${meta.id}/*.genotype.json")     , emit: genotype
    // tuple val(meta), path("${meta.id}/*.genotype.log")      , emit: genotype_log
    // tuple val(meta), path("${meta.id}/*.alignment.p")       , emit: aligment

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "0.5.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mv $fq /usr/local/share/arcas-hla-0.5.0-1
    
    cd /usr/local/share/arcas-hla-0.5.0-1

    arcasHLA reference --version 3.24.0

    arcasHLA \
        genotype \
        $args \
        -t $task.cpus \
        -o ${prefix} \
        $fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arcashla: $VERSION
    END_VERSIONS
    """
}
