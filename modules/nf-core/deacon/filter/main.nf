process DEACON_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deacon:0.10.0--h4349ce8_0':
        'biocontainers/deacon:0.10.0--h4349ce8_0' }"


    input:
    tuple val(meta), path(index)
    tuple val(meta2), path(reads)
    val (save_summary)

    output:
    tuple val(meta2), path("filtered_reads/*.fastq.gz"), emit: filtered_reads
    tuple val(meta2), path("*.json")                   , emit: summary_json, optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_input = meta2.single_end ? reads : "${reads.sort()[0]} ${reads.sort()[1]}"
    def summary_arg = save_summary ? "-s ${prefix}_summary.json" : ""
    def output_arg = meta2.single_end ? "-o filtered_reads/${prefix}.fastq.gz" : "-o filtered_reads/${prefix}_1.fastq.gz -O filtered_reads/${prefix}_2.fastq.gz"
    """
    mkdir -p filtered_reads/

    deacon \\
        filter \\
        $args \\
        --threads ${task.cpus} \\
        $index \\
        $reads_input \\
        $summary_arg \\
        $output_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    """
    mkdir -p filtered_reads
    if [[ "${meta2.single_end}" == "true" ]]; then
        echo "" | gzip -c > filtered_reads/${prefix}.fastq.gz
    else
        echo "" | gzip -c > filtered_reads/ ${prefix}_1.fastq.gz
        echo "" | gzip -c > filtered_reads/ ${prefix}_2.fastq.gz
    fi

    if [[ "${save_summary}" == "true" ]]; then
        touch ${prefix}_summary.json
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deacon: \$(deacon --version | head -n1 | sed 's/deacon //g')
    END_VERSIONS
    """
}
