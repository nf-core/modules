
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process MCSTAGING_PHENOIMAGER2MC {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/schapirolabor/phenoimager2mc:v0.1.1" // how to specify a version?

    input:
    tuple val(meta) , path(folder)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.tif"), emit: tif

    when:
    task.ext.when == null || task.ext.when

    script:
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Phenoimager2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 /phenoimager2mc/scripts/phenoimager2mc.py \
        -i ${folder} \
        -o "${prefix}.tif" \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phenoimager2mc: \$(python3 /phenoimager2mc/scripts/phenoimager2mc.py --version | sed 's/v//g')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Phenoimager2mc module in conda does not exist. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch input
    touch "${prefix}.tif"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phenoimager2mc: \$(python3 /scripts/phenoimager2mc.py --version | sed 's/v//g')
    END_VERSIONS
    """
}
