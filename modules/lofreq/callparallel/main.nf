process LOFREQ_CALLPARALLEL {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::lofreq=2.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4"
    } else {
        container "quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    lofreq \\
        call-parallel \\
        --pp-threads $task.cpus \\
        $options.args \\
        -f $fasta \\
        -o ${prefix}.vcf.gz \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}
