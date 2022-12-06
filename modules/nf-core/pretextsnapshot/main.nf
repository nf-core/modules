process PRETEXTSNAPSHOT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::pretextsnapshot=0.0.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextsnapshot:0.0.4--h7d875b9_0':
        'quay.io/biocontainers/pretextsnapshot:0.0.4--h7d875b9_0' }"

    input:
    tuple val(meta), path(pretext_map)

    output:
    tuple val(meta), path('*.{jpg,png,bmp}', glob: true), emit: image
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = task.ext.outdir ?: '.'
    """
    PretextSnapshot \\
        $args \\
        -m $pretext_map \\
        --prefix $prefix \\
        -o $outdir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextsnapshot: \$(echo \$(PretextSnapshot --version 2>&1) | sed 's/^.*PretextSnapshot Version //' )
    END_VERSIONS
    """
}
