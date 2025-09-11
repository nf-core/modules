process HOSTILE_CLEAN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hostile:2.0.2--pyhdfd78af_0'
        : 'biocontainers/hostile:2.0.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(reference_name), path(reference_dir)

    output:
    tuple val(meta), path('cleaned_reads/*.fastq.gz'), emit: fastq
    tuple val(meta), path('*.json'), emit: json
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_cmd = meta.single_end ? "--fastq1 ${reads.sort()[0]}" : "--fastq1 ${reads.sort()[0]} --fastq2 ${reads.sort()[1]}"
    """
    export HOSTILE_CACHE_DIR=${reference_dir}
    mkdir cleaned_reads/

    ## Reorder the reads for reproducibility
    ## Set offline as we never want this process to auto-download reference files as required input channel
    hostile \\
        clean \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index ${reference_dir}/${reference_name} \\
        --output cleaned_reads/ \\
        --reorder \\
        --airplane \\
        | tee > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_cmd = meta.single_end ? "echo '' | gzip -c > cleaned_reads/${prefix}.clean_2.fastq.gz" : ""

    """
    echo "hostile \\
        clean \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index ${reference_dir}/${reference_name} \\
        --output cleaned_reads/ \\
        --reorder \\
        --airplane \\
        | tee > ${prefix}.json"

    export HOSTILE_CACHE_DIR=${reference_dir}
    mkdir cleaned_reads/
    echo "" | gzip -c > cleaned_reads/${prefix}.clean_1.fastq.gz
    ${reads_cmd}
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """
}
