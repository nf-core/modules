process GLIMPSE2_SPLITREFERENCE {
    tag "${meta.id}"
    label 'process_low'

    beforeScript """
    if cat /proc/cpuinfo | grep avx2 -q
    then
        echo "Feature AVX2 present"
    else
        echo "Feature AVX2 not present on node"
        exit 1
    fi
    """

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.1--h46b9e50_1'
        : 'biocontainers/glimpse-bio:2.0.1--h46b9e50_1'}"

    input:
    tuple val(meta), path(reference), path(reference_index), val(input_region), val(output_region), path(map)

    output:
    tuple val(meta), path("*.bin"), emit: bin_ref
    tuple val("${task.process}"), val('glimpse2'), eval("GLIMPSE2_split_reference --help | grep -oE 'v[0-9.]+' | cut -c2-"), topic: versions, emit: versions_glimpse2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}_${output_region.replace(":", "_")}"
    def map_command = map             ? "--map ${map}" : ""

    """
    GLIMPSE2_split_reference \\
        ${args} \\
        --reference ${reference} \\
        ${map_command} \\
        --input-region ${input_region} \\
        --output-region ${output_region} \\
        --thread ${task.cpus} \\
        --output ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${output_region.replace(":", "_")}"
    """
    touch ${prefix}.bin
    """
}
