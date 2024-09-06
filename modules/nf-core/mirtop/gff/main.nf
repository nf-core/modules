process MIRTOP_GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' :
        'biocontainers/mulled-v2-0c13ef770dd7cc5c76c2ce23ba6669234cf03385:63be019f50581cc5dfe4fc0f73ae50f2d4d661f7-0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(hairpin)
    tuple val(meta3), path(gtf), val(species)

    output:
    tuple val(meta), path("mirtop/${bam.baseName}.gff")  , emit: sample_gff
    tuple val(meta), path("mirtop/mirtop.gff")           , emit: mirtop_gff
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mirtop \\
        gff \\
        $args \\
        --sps $species \\
        --hairpin $hairpin \\
        --gtf $gtf \\
        -o mirtop \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir mirtop
    touch mirtop/mirtop.gff
    touch mirtop/sim_isomir_sort.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """
}
