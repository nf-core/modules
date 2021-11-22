process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0'
    } else {
        container 'quay.io/biocontainers/bowtie2:2.4.4--py36hd4290be_0'
    }

    input:
    path fasta

    output:
    path 'bowtie2'      , emit: index
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build $args --threads $task.cpus $fasta bowtie2/${fasta.baseName}
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
