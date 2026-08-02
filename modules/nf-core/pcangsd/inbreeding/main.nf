process PCANGSD_INBREEDING {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pcangsd%3A1.36.4--py313h5d164f8_1':
        'quay.io/biocontainers/pcangsd:1.36.4--py313h5d164f8_1' }"

    input:
    tuple val(meta), path(beagle_file)

    output:
    tuple val(meta), path("*.inbreed.samples"), emit: inbreeding_coefficients
    tuple val("${task.process}"), val('pcangsd'), eval("pcangsd --version"), topic: versions, emit: versions_pcangsd

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pcangsd \\
        --threads ${task.cpus} \\
        --beagle ${beagle_file} \\
        --inbreed-samples \\
        --out ${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.inbreed.samples
    """
}
