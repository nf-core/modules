process UNZIP {
    tag "$archive"
    label 'process_single'

    conda "conda-forge::p7zip=16.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:15.09--h2d50403_4' :
        'quay.io/biocontainers/p7zip:15.09--h2d50403_4' }"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}/")        , emit: unzipped_archive
    tuple val(meta), path("${prefix}_files/**"), emit: files, optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    if ( archive instanceof List && archive.name.size > 1 ) { exit 1, "[UNZIP] error: 7za only accepts a single archive as input. Please check module input." }
    
    prefix = task.ext.prefix ?: ( meta.id ? "${meta.id}" : archive.baseName)
    def emit_files = args2.contains("--emit-files") ? "true": ''
    """
    7za \\
        x \\
        -o"${prefix}"/ \\
        $args \\
        $archive
    
    if [[ -n "$emit_files" ]]; then
        ln -s $prefix ${prefix}_files
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
