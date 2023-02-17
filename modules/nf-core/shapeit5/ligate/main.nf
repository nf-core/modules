process SHAPEIT5_LIGATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::shapeit5=1.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shapeit5:1.0.0--h0c8ee15_0':
        'quay.io/biocontainers/shapeit5:1.0.0--h0c8ee15_0'}"

    input:
    tuple val(meta), path(input_list), path (input_list_index)

    output:
    tuple val(meta), path("*.{vcf,bcf,vcf.gz,bcf.gz}"), emit: merged_variants
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "vcf.gz"
    """
    printf "%s\\n" $input_list | tr -d '[],' > all_files.txt

    SHAPEIT5_ligate \\
        $args \\
        --input all_files.txt \\
        --thread $task.cpus \\
        --output ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            shapeit5: "\$(SHAPEIT5_ligate | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -n 1)"
    END_VERSIONS
    """
}
