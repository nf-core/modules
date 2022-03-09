
process ARIMA_FILTERFIVEEND {
    tag "$meta.id"
    label 'process_medium'
    def version = '0.001-c24'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the arima filter_five_end_v1 tool. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/arima:${version}"

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        samtools view -h ${read} | filter_five_end.pl | samtools view -Sb $args -@ $task.cpus -o ${prefix}_filtered.bam -
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
            arima: ${version}
        END_VERSIONS
        """
    } else {
        """
        samtools view -h ${read[0]} | filter_five_end.pl | samtools view -Sb $args -@ $task.cpus -o ${prefix}_1_filtered.bam -
        samtools view -h ${read[1]} | filter_five_end.pl | samtools view -Sb $args -@ $task.cpus -o ${prefix}_2_filtered.bam -
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
            arima: ${version}
        END_VERSIONS
        """
    }
}
