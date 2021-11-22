process MALTEXTRACT {

    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hops=0.35" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hops:0.35--hdfd78af_1' :
        'quay.io/biocontainers/hops:0.35--hdfd78af_1' }"

    input:
    path rma6
    path taxon_list
    path ncbi_dir

    output:
    path "results"      , emit: results
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    MaltExtract \\
        -Xmx${task.memory.toGiga()}g \\
        -p $task.cpus \\
        -i ${rma6.join(' ')} \\
        -t $taxon_list \\
        -r $ncbi_dir \\
        -o results/ \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(MaltExtract --help | head -n 2 | tail -n 1 | sed 's/MaltExtract version//')
    END_VERSIONS
    """
}
