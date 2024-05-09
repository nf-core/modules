process SCSPLIT_COUNT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'irenerobles93/scsplit:1.0.8'

    input:
    tuple val(meta), path(bam), path(bai), path(vcf), path(barcode)

    output:
    tuple val(meta), path('*_alt_filtered.csv'), emit: alt_filtered
    tuple val(meta), path('*_ref_filtered.csv'), emit: ref_filtered
    tuple val(meta), path('*_scSplit.log')     , emit: log
    path 'versions.yml'                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scSplit count -v $vcf -i $bam -b $barcode -r ${prefix}_ref_filtered.csv -a ${prefix}_alt_filtered.csv -o .
    mv scSplit.log ${prefix}_scSplit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scsplit: 1.0.8
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_ref_filtered.csv
    touch ${prefix}_alt_filtered.csv
    touch ${prefix}_scSplit.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scsplit: 1.0.8
    END_VERSIONS
    """
}
