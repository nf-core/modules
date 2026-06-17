process SHAPEIT5_SWITCH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fa/fa06870e893f9045944461e5b674adffec6deec41e8496b63ec54d44a98d6134/data'
        : 'community.wave.seqera.io/library/shapeit5:5.1.1--09a6cb254ece8f6e'}"

    input:
    tuple val(meta), path(estimate), path(estimate_index), val(region), path(pedigree), path(truth), path(truth_index), path(freq), path(freq_index)

    output:
    tuple val(meta), path("*.txt.gz"), emit: errors
    tuple val("${task.process}"), val('shapeit5'), eval('SHAPEIT5_switch | sed "5!d;s/^.*Version *: //; s/ .*$//"'), topic: versions, emit: versions_shapeit5


    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def freq_cmd     = freq     ? "--frequency ${freq}"    : ""
    def pedigree_cmd = pedigree ? "--pedigree ${pedigree}" : ""

    """
    SHAPEIT5_switch \\
        ${args} \\
        --estimation ${estimate} \\
        --region ${region} \\
        --validation ${truth} \\
        ${freq_cmd} \\
        ${pedigree_cmd} \\
        --thread ${task.cpus} \\
        --output ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.block.switch.txt.gz
    echo "" | gzip > ${prefix}.calibration.switch.txt.gz
    echo "" | gzip > ${prefix}.flipsAndSwitches.txt.gz
    echo "" | gzip > ${prefix}.frequency.switch.txt.gz
    echo "" | gzip > ${prefix}.sample.switch.txt.gz
    echo "" | gzip > ${prefix}.sample.typing.txt.gz
    echo "" | gzip > ${prefix}.type.switch.txt.gz
    echo "" | gzip > ${prefix}.variant.switch.txt.gz
    echo "" | gzip > ${prefix}.variant.typing.txt.gz
    """
}
