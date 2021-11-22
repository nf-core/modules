process SAMTOOLS_BAM2FQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(meta), path(inputbam)
    val(split)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (split){
        """
        samtools \\
            bam2fq \\
            $options.args \\
            -@ $task.cpus \\
            -1 ${prefix}_1.fq.gz \\
            -2 ${prefix}_2.fq.gz \\
            -0 ${prefix}_other.fq.gz \\
            -s ${prefix}_singleton.fq.gz \\
            $inputbam

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        samtools \\
            bam2fq \\
            $options.args \\
            -@ $task.cpus \\
            $inputbam >${prefix}_interleaved.fq.gz

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}
