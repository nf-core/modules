process SHAPEIT5_PHASERARE {
    tag "$meta.id"
    label 'process_low'

    beforeScript  """
    if cat /proc/cpuinfo | grep avx2 -q
    then
        echo "Feature AVX2 present on host"
    else
        echo "Feature AVX2 not present on host"
        exit 1
    fi
    """

    conda "bioconda::shapeit5=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:1.0.0--h0c8ee15_0':
        'biocontainers/shapeit5:1.0.0--h0c8ee15_0' }"

    input:
        tuple val(meta) , path(input_plain), path(input_plain_index), path(pedigree), val(input_region)
        tuple val(meta2), path(scaffold)   , path(scaffold_index)   , val(scaffold_region)
        tuple val(meta3), path(map)

    output:
        tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: phased_variant
        path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    if ("$input_plain" == "${prefix}.${suffix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def map_command       = map       ? "--map $map"             : ""
    def pedigree_command  = pedigree  ? "--pedigree $pedigree"   : ""

    """
    SHAPEIT5_phase_rare \\
        $args \\
        --input-plain $input_plain \\
        --scaffold $scaffold \\
        $map_command \\
        $pedigree_command \\
        --input-region $input_region \\
        --scaffold-region $scaffold_region \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shapeit5: "\$(SHAPEIT5_phase_rare | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
