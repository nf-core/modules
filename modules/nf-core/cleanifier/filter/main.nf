process CLEANIFIER_FILTER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cleanifier:1.3.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/cleanifier:1.3.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index), path(info)

    output:
    tuple val(meta), path("*.clean*keep.fastq.gz"), emit: reads
    tuple val(meta), path("*.clean*filter.fastq.gz"), emit: removed, optional: true
    tuple val("${task.process}"), val('cleanifier'), eval("cleanifier --version"), emit: versions_cleanifier, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--fastq ${reads}" : "--fastq ${reads[0]} --pairs ${reads[1]}"
    def index_base = index.baseName
    def threads_compression = Math.max(1, (task.cpus / 4) as int)
    def threads = Math.max(1, (task.cpus - threads_compression) as int)
    """
    cleanifier filter \\
        ${input_reads} \\
        --index ${index_base} \\
        --out ${prefix}.clean \\
        --threads ${threads} \\
        --compression gz \\
        --compression-threads ${threads_compression} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def keep_host = args.contains("--keep-host")
    if (meta.single_end) {
        """
        echo "" | gzip > ${prefix}.clean_keep.fastq.gz
        ${keep_host ? "echo \"\" | gzip > ${prefix}.clean_filter.fastq.gz" : ""}
        """
    }
    else {
        """
        echo "" | gzip > ${prefix}.clean_1_keep.fastq.gz
        echo "" | gzip > ${prefix}.clean_2_keep.fastq.gz
        ${keep_host ? "echo \"\" | gzip > ${prefix}.clean_1_filter.fastq.gz" : ""}
        ${keep_host ? "echo \"\" | gzip > ${prefix}.clean_2_filter.fastq.gz" : ""}
        """
    }
}
