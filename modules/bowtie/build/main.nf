process BOWTIE_BUILD {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hed8969a_1'
    } else {
        container 'quay.io/biocontainers/bowtie:1.3.0--py38hed8969a_1'
    }

    input:
    path fasta

    output:
    path 'bowtie'       , emit: index
    path "versions.yml" , emit: versions

    script:
    """
    mkdir bowtie
    bowtie-build --threads $task.cpus $fasta bowtie/${fasta.baseName}
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
