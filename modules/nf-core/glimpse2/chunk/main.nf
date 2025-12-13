process GLIMPSE2_CHUNK {
    tag "${meta.id}"
    label 'process_low'

    beforeScript """
        if cat /proc/cpuinfo | grep avx2 -q
        then
            echo "Feature AVX2 present on host"
        else
            echo "Feature AVX2 not present on host"
            exit 1
        fi
    """
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.1--h46b9e50_1'
        : 'biocontainers/glimpse-bio:2.0.1--h46b9e50_1'}"

    input:
    tuple val(meta), path(input), path(input_index), val(region), path(map)
    val model

    output:
    tuple val(meta), path("*.txt"), emit: chunk_chr
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ""
    def map_cmd = map ? "--map ${map}" : ""

    """
    GLIMPSE2_chunk \\
        ${args} \\
        ${map_cmd} \\
        --${model} \\
        --input ${input} \\
        --region ${region} \\
        --threads ${task.cpus} \\
        --output ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_chunk --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "${meta.id}\t${region}\t0\t0\t0\t0\t0\t0" > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse2: "\$(GLIMPSE2_chunk --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
