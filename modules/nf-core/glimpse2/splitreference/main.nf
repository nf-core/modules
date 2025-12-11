process GLIMPSE2_SPLITREFERENCE {
    tag "$meta.id"
    label 'process_low'

    beforeScript  """
    if cat /proc/cpuinfo | grep avx2 -q
    then
        echo "Feature AVX2 present"
    else
        echo "Feature AVX2 not present on node"
        exit 1
    fi
    """

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.1--h46b9e50_1':
        'biocontainers/glimpse-bio:2.0.1--h46b9e50_1' }"

    input:
        tuple val(meta) , path(reference), path(reference_index), val(input_region), val(output_region), path(map)


    output:
        tuple val(meta), path("*.bin"), emit: bin_ref
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}_${output_region.replace(":","_")}"
    def map_command = map             ? "--map $map" : ""

    """
    GLIMPSE2_split_reference \\
        $args \\
        --reference $reference \\
        $map_command \\
        --input-region $input_region \\
        --output-region $output_region \\
        --thread $task.cpus \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_split_reference --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${output_region.replace(":","_")}"
    """
    touch ${prefix}.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_split_reference --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
