process METABINNER_GENCOVERAGEFILE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8dab51352a1e2701daf87620104c5c6e7148811ee132769f299e8fa434fc194e/data'
:         'community.wave.seqera.io/library/metabinner:1.4.4--60cf94125f522897' }"

    input:
    tuple val(meta), path(assembly), path(reads, stageAs: 'reads/*')

    output:
    tuple val(meta), path("${prefix}_coverage_profile.tsv"), emit: coverage
    tuple val(meta), path("${prefix}_coverage_profile_f1k.tsv"), emit: coverage_f1k
    tuple val(meta), path("${prefix}_depth.txt"), emit: depth
    tuple val(meta), path("${prefix}_s*.bam"), emit: bam
    // hardcoded: gen_coverage_file.sh has no --version flag
    tuple val("${task.process}"), val('metabinner'), val('1.4.4'), topic: versions, emit: versions_metabinner

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads_list = [reads].flatten()
    def interleaved = args.contains('--interleaved')
    def paired = !meta.single_end && !interleaved
    def read_mode = paired ? '-f _1.fastq.gz -r _2.fastq.gz' : (interleaved ? '' : '--single-end')
    def link_reads = paired
        ? reads_list.collate(2).withIndex().collect { pair, i -> "ln -s ${pair[0]} ${prefix}_s${i + 1}_1.fastq.gz; ln -s ${pair[1]} ${prefix}_s${i + 1}_2.fastq.gz" }.join('\n    ')
        : reads_list.withIndex().collect { file, i -> "ln -s ${file} ${prefix}_s${i + 1}.fastq" }.join('\n    ')
    def read_files = paired
        ? reads_list.collate(2).indices.collect { i -> "${prefix}_s${i + 1}_1.fastq.gz ${prefix}_s${i + 1}_2.fastq.gz" }.join(' ')
        : reads_list.indices.collect { i -> "${prefix}_s${i + 1}.fastq" }.join(' ')
    """
    ${link_reads}

    gzip -cdf ${assembly} > ${prefix}_assembly.fa

    gen_coverage_file.sh \\
        -t ${task.cpus} \\
        -a ${prefix}_assembly.fa \\
        -o ${prefix}_out \\
        -b . \\
        ${read_mode} \\
        ${args} \\
        ${read_files}

    mv ${prefix}_out/coverage_profile.tsv ${prefix}_coverage_profile.tsv
    mv ${prefix}_out/coverage_profile_f1k.tsv ${prefix}_coverage_profile_f1k.tsv
    mv ${prefix}_out/work_files/mb2_master_depth.txt ${prefix}_depth.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"

    touch ${prefix}_coverage_profile.tsv
    touch ${prefix}_coverage_profile_f1k.tsv
    touch ${prefix}_depth.txt
    touch ${prefix}_s1.bam
    """
}
