
process AMPLIFY_PREDICT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::amplify=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/amplify:1.0.3--py36hdfd78af_0':
        'quay.io/biocontainers/amplify:1.0.3--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(faa)
    
    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml",            emit: versions

    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    amplify\\
     -s "${faa}"\\
     $args
    mv *.tsv amplify_output.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amplify_predict: 1.0.3
    END_VERSIONS
    """
}
