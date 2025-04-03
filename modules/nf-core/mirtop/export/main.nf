process MIRTOP_EXPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/mirtop_pybedtools_pysam_samtools_pruned:60b8208f3dbb2910"

    input:
    tuple val(meta), path(mirtop_gff)
    tuple val(meta2), path(hairpin)
    tuple val(meta3), path(gtf), val(species)

    output:
    tuple val(meta), path("export/*_rawData.tsv")     , emit: tsv, optional: true
    tuple val(meta), path("export/*.fasta")           , emit: fasta, optional: true
    tuple val(meta), path("export/*.vcf*")            , emit: vcf  , optional: true
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--format isomir'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mirtop \\
        export \\
        $args \\
        --hairpin $hairpin\\
        --gtf $gtf \\
        --sps $species \\
        -o export \\
        $mirtop_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir export
    touch export/${prefix}.fasta
    touch export/${prefix}.vcf
    touch export/${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mirtop: \$(echo \$(mirtop --version 2>&1) | sed 's/^.*mirtop //')
    END_VERSIONS
    """
}
