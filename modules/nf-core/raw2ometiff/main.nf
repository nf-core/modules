
process RAW2OMETIFF {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/blosc_raw2ometiff:6ae5a6e69f313aa5"

    input:
    tuple val(meta), path(zarr_dir)

    output:
    tuple val(meta), path("*.tiff"), emit: tiff
    path "versions.yml"           ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    raw2ometiff \\
        ${zarr_dir} \\
        ${prefix}.tiff \\
        --max_workers $task.cpus \\
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raw2ometiff: \$(raw2ometiff --version |& sed -n '1s/Version = //p')
        bio-formats: \$(bioformats2raw --version |& sed -n '2s/Bio-Formats version = //p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${prefix}.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raw2ometiff: \$(raw2ometiff --version |& sed -n '1s/Version = //p')
        bio-formats: \$(bioformats2raw --version |& sed -n '2s/Bio-Formats version = //p')
    END_VERSIONS
    """
}
