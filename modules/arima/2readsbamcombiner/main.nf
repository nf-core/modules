process ARIMA_2READSBAMCOMBINER {
    tag "$meta.id"
    label 'process_low'

    def version = '0.001-c24'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the arima two reads bam combiner. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/arima:${version}"

    input:
    tuple val(meta), path(bams)
    val qscore

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    two_read_bam_combiner.pl $bams 'samtools' $qscore | samtools sort $args -@ $task.cpus -o ${prefix}_combined.bam -T $prefix -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        arima: $version
    END_VERSIONS
    """
}
