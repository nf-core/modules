process HOSTILE_CLEAN {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7caca3a47606de8e3460b35823193a471272aa6ab7cfafbf9aabf4615c9fa181/data'
        : 'community.wave.seqera.io/library/hostile:2.0.2--a7f5e5d341b6b94b'}"

    input:
    tuple val(meta)          , path(reads, stageAs: "input_reads/")
    tuple val(reference_name), path(reference_dir)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    tuple val(meta), path('*.json')    , emit: json
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { it.simpleName }
    def reads_cmd    = meta.single_end ? "--fastq1 ${sorted_reads[0]}" : "--fastq1 ${sorted_reads[0]} --fastq2 ${sorted_reads[1]}"
    """
    export HOSTILE_CACHE_DIR=${reference_dir}

    ## Reorder the reads for reproducibility
    ## Set offline as we never want this process to auto-download reference files as required input channel
    hostile \\
        clean \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index ${reference_name} \\
        --output . \\
        --reorder \\
        --airplane \\
        | tee > ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """

    stub:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { it.simpleName }
    def reads_cmd    = meta.single_end ? "--fastq1 ${sorted_reads[0]}" : "--fastq1 ${sorted_reads[0]} --fastq2 ${sorted_reads[1]}"
    def fake_read2   = !meta.single_end ? "echo '' | gzip -c > ${prefix}.clean_2.fastq.gz" : ""
    """
    export HOSTILE_CACHE_DIR=${reference_dir}

    echo "hostile \\
        clean \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads_cmd} \\
        --index ${reference_name} \\
        --output . \\
        --reorder \\
        --airplane \\
        | tee > ${prefix}.json"

    export HOSTILE_CACHE_DIR=${reference_dir}
    echo "" | gzip -c > ${prefix}.clean_1.fastq.gz
    ${fake_read2}

    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version)
    END_VERSIONS
    """
}
