process ALEVINFRY_GENERATEPERMITLIST {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::alevin-fry=0.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/alevin-fry%3A0.5.0--h9f5acd7_1':
        'quay.io/biocontainers/alevin-fry:0.5.0--h9f5acd7_1' }"

    input:
    tuple val(meta), path(inputdir)

    output:
    path val(meta), path("${prefix}_af_permit_list"), emit: permitlist 
    tuple val(meta), path("${prefix}_af_permit_list/generate_permit_list.json") , emit: log
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    alevin-fry generate-permit-list \\
        --input $inputdir \\
        $args \\
        --output-dir ${prefix}_af_permit_list

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevinfry: \$(echo \$(alevin-fry --version 2>&1) | sed 's/^alevin-fry //g' )
    END_VERSIONS
    """
}