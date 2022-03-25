process METAMAPS_CLASSIFY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::metamaps=0.1.98102e9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metamaps:0.1.98102e9--h176a8bc_0':
        'quay.io/biocontainers/metamaps:0.1.98102e9--h176a8bc_0' }"

    input:
    tuple val(meta), path(classification_res)
    path database_folder

    output:
    tuple val(meta), path("*classification_res")                           , emit: classification_res
    tuple val(meta), path("*classification_res.meta")                      , emit: meta_file
    tuple val(meta), path("*classification_res.meta.unmappedReadsLengths") , emit: meta_unmappedReadsLengths
    tuple val(meta), path("*classification_res.parameters")                , emit: para_file

    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    metamaps \\
        classify \\
        --mappings $classification_res \\
        --threads $task.cpus \\
        --DB $database_folder \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metamaps: \$(echo \$(metamaps 2>&1) | grep -E "^MetaMaps\sv\s" | sed -E 's/^MetaMaps\sv\s([0-9]+\.[0-9]+)/\1/')
    END_VERSIONS
    """
}
