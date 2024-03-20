process DATAVZRD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/datavzrd:2.23.2':
        'biocontainers/datavzrd:2.23.2' }"

    input:
    path config
    path output_path


    output:
    path "*.html", emit: output_path
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def overwrite = overwrite ? "--overwrite-output" : ""
    def debug = debug ? "--debug" : ""
    def webview_url = webview_url ? "--webview-url ${}" : ""

    """
    datavzrd \\
        ${debug} \\
        ${overwrite} \\
        ${debug} \\
        ${webview_url} \\
        --output ${input.output_path} \\
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd version | sed -e "s/datavzrd v//g" ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def overwrite = overwrite ? "--overwrite-output" : ""
    def debug = debug ? "--debug" : ""
    def webview_url = webview_url ? "--webview-url ${}" : ""
    
    """
    datavzrd \\
        ${debug} \\
        ${overwrite} \\
        ${debug} \\
        ${webview_url} \\
        --output ${input.output_path} \\
      

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        datavzrd: \$(echo \$( datavzrd version | sed -e "s/datavzrd v//g" ))
    END_VERSIONS
    """
}
