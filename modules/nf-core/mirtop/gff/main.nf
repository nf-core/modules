process MIRTOP_GFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0da43138fd5dfa0d365ef64ba39061102efa11256aea303791869ce46044a3df/data':
        'community.wave.seqera.io/library/mirtop_pybedtools_pysam_samtools:b9705c2683e775b8' }"

    input:
    tuple val(meta), path(bam, arity:'1..*')
    tuple val(meta2), path(hairpin)
    tuple val(meta3), path(gtf), val(species)

    output:
    tuple val(meta), path("mirtop/*mirtop.gff")           , emit: gff
    path "versions.yml"                                   , emit: versions

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

    mv mirtop/mirtop.gff mirtop/${prefix}_mirtop.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir mirtop
    touch mirtop/${prefix}_mirtop.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """
}
