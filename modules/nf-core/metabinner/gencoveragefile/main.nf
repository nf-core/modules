process METABINNER_GENCOVERAGEFILE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8dab51352a1e2701daf87620104c5c6e7148811ee132769f299e8fa434fc194e/data'
        : 'community.wave.seqera.io/library/metabinner:1.4.4--60cf94125f522897'}"

    input:
    tuple val(meta), path(assembly), path(reads, stageAs: 'reads/?/*')

    output:
    tuple val(meta), path("${prefix}_coverage_profile.tsv"), emit: coverage
    tuple val(meta), path("${prefix}_coverage_profile_f1k.tsv"), emit: coverage_f1k
    tuple val(meta), path("${prefix}_depth.txt"), emit: depth
    tuple val(meta), path("${prefix}_s*.bam"), emit: bam
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('metabinner'), val('1.4.4'), topic: versions, emit: versions_metabinner

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads_list = [reads].flatten()
    def interleaved = args.contains('--interleaved')
    if (args.contains('--single-end')) {
        error("METABINNER_GENCOVERAGEFILE: set meta.single_end instead of passing --single-end in ext.args")
    }
    if (!interleaved && meta.single_end == null) {
        error("METABINNER_GENCOVERAGEFILE: meta.single_end must be set unless --interleaved is passed in ext.args")
    }
    def paired = !meta.single_end && !interleaved
    if (paired && reads_list.size() % 2 != 0) {
        error("METABINNER_GENCOVERAGEFILE: paired-end mode expects an even number of read files, got ${reads_list.size()}")
    }
    def samples = reads_list.collate(paired ? 2 : 1)
    // Zero-padded so that the lexically sorted *.bam glob, and hence the coverage columns, follow input order
    def width = Math.max(3, samples.size().toString().length())
    def names = samples.indices.collect { i -> "${prefix}_s${String.format("%0${width}d", i + 1)}" }
    def link_reads = [samples, names].transpose().collect { files, name ->
        paired
            ? "ln -s ${files[0]} ${name}_1.fastq.gz; ln -s ${files[1]} ${name}_2.fastq.gz"
            : "ln -s ${files[0]} ${name}.fastq"
    }.join('\n    ')
    def read_files = names.collect { name -> paired ? "${name}_1.fastq.gz ${name}_2.fastq.gz" : "${name}.fastq" }.join(' ')
    def read_mode = paired ? '-f _1.fastq.gz -r _2.fastq.gz' : (interleaved ? '' : '--single-end')
    // args go first: the upstream parser keeps the last value of a repeated option
    """
    ${link_reads}

    gen_coverage_file.sh \\
        ${args} \\
        -t ${task.cpus} \\
        -a ${assembly} \\
        -o ${prefix}_out \\
        -b . \\
        ${read_mode} \\
        ${read_files}

    mv ${prefix}_out/coverage_profile.tsv ${prefix}_coverage_profile.tsv
    # 1.4.4 ignores -l and always filters at 1000 bp; later releases write coverage_profile_f<len>.tsv instead
    mv ${prefix}_out/coverage_profile_f1k.tsv ${prefix}_coverage_profile_f1k.tsv
    mv ${prefix}_out/work_files/mb2_master_depth.txt ${prefix}_depth.txt
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads_list = [reads].flatten()
    def paired = !meta.single_end && !args.contains('--interleaved')
    def samples = reads_list.collate(paired ? 2 : 1)
    def width = Math.max(3, samples.size().toString().length())
    def touch_bams = samples.indices.collect { i -> "touch ${prefix}_s${String.format("%0${width}d", i + 1)}.bam" }.join('\n    ')
    """
    echo "${args}"

    touch ${prefix}_coverage_profile.tsv
    touch ${prefix}_coverage_profile_f1k.tsv
    touch ${prefix}_depth.txt
    ${touch_bams}
    """
}
