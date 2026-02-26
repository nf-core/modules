process SHAPEIT5_PHASERARE {
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
        ? 'https://depot.galaxyproject.org/singularity/shapeit5:5.1.1--hb60d31d_0'
        : 'biocontainers/shapeit5:5.1.1--hb60d31d_0'}"

    input:
    tuple val(meta), path(input), path(input_index), path(pedigree), val(input_region), path(scaffold), path(scaffold_index), val(scaffold_region), path(map)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: phased_variant
    tuple val("${task.process}"), val('shapeit5'), eval('SHAPEIT5_phase_rare | sed "5!d;s/^.*Version *: //; s/ .*$//"'), topic: versions, emit: versions_shapeit5

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    if ("${input}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    if ("${scaffold}" == "${prefix}.${suffix}") {
        error("Scaffold and output names are the same, set prefix in module configuration to disambiguate!")
    }

    def map_command      = map      ? "--map ${map}"           : ""
    def pedigree_command = pedigree ? "--pedigree ${pedigree}" : ""

    """
    SHAPEIT5_phase_rare \\
        ${args} \\
        --input ${input} \\
        --scaffold ${scaffold} \\
        ${map_command} \\
        ${pedigree_command} \\
        --input-region ${input_region} \\
        --scaffold-region ${scaffold_region} \\
        --thread ${task.cpus} \\
        --output ${prefix}.${suffix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    def create_cmd = suffix.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    """
    ${create_cmd} ${prefix}.${suffix}
    """
}
