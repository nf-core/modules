process PLASMIDID {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::plasmidid=1.6.5' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasmidid:1.6.5--hdfd78af_0' :
        'quay.io/biocontainers/plasmidid:1.6.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(scaffold)
    path  fasta

    output:
    tuple val(meta), path("${prefix}/*final_results.html"), emit: html
    tuple val(meta), path("${prefix}/*final_results.tab") , emit: tab
    tuple val(meta), path("${prefix}/images/")            , emit: images
    tuple val(meta), path("${prefix}/logs/")              , emit: logs
    tuple val(meta), path("${prefix}/data/")              , emit: data
    tuple val(meta), path("${prefix}/database/")          , emit: database
    tuple val(meta), path("${prefix}/fasta_files/")       , emit: fasta_files
    tuple val(meta), path("${prefix}/kmer/")              , emit: kmer
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    plasmidID \\
        -d $fasta \\
        -s $prefix \\
        -c $scaffold \\
        $args \\
        -o .

    mv NO_GROUP/$prefix ./$prefix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmidid: \$(echo \$(plasmidID --version 2>&1))
    END_VERSIONS
    """
}
