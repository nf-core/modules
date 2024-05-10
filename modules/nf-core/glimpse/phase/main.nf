process GLIMPSE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--hce55b13_1':
        'biocontainers/glimpse-bio:1.1.1--hce55b13_1' }"

    input:
        tuple val(meta) , path(input), path(input_index), path(samples_file), val(input_region), val(output_region), path(reference), path(reference_index), path(map)

    output:
        tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: phased_variants
        path "versions.yml"                               , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${input_region.replace(":","_")}"
    def suffix = task.ext.suffix ?: "vcf.gz"

    def map_command           = map           ? "--map $map"                    :""
    def samples_file_command  = samples_file  ? "--samples-file $samples_file"  :""

    """
    GLIMPSE_phase \\
        $args \\
        --input $input \\
        --reference $reference \\
        $map_command \\
        $samples_file_command \\
        --input-region $input_region \\
        --output-region $output_region \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: "\$(GLIMPSE_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${input_region.replace(":","_")}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: "\$(GLIMPSE_phase --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """
}
