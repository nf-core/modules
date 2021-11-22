process SEQTK_SAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::seqtk=1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
    } else {
        container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
    }

    input:
    tuple val(meta), path(reads)
    val sample_size

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        seqtk \\
            sample \\
            $args \\
            $reads \\
            $sample_size \\
            | gzip --no-name > ${prefix}.fastq.gz \\

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        if (!(options.args ==~ /.*-s[0-9]+.*/)) {
            options.args = options.args + " -s100"
        }
        """
        seqtk \\
            sample \\
            $args \\
            ${reads[0]} \\
            $sample_size \\
            | gzip --no-name > ${prefix}_1.fastq.gz \\

        seqtk \\
            sample \\
            $args \\
            ${reads[1]} \\
            $sample_size \\
            | gzip --no-name > ${prefix}_2.fastq.gz \\

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
