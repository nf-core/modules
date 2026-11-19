process SAMTOOLS_MULTICOMMAND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c2f02ceb8cbe5aa37930dcbc0de6f1d07dfb292178c7e780c44ab2f556bf829f/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--84fe6b04a60dfab2'}"

    input:
    tuple val(meta), path(input, arity: '1..*'), path(index, arity: '0..*')
    tuple val(meta2), path(fasta), path(fai)
    val chain

    output:
    // Alignment format outputs (view, sort, markdup, fixmate, merge, cat, collate)
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.cram"), optional: true, emit: cram
    tuple val(meta), path("*.sam"), optional: true, emit: sam
    tuple val(meta), path("*.{bai,csi,crai}"), optional: true, emit: index

    // Sequence outputs (fasta, fastq)
    tuple val(meta), path("*_{1,2}.fasta.gz"), optional: true, emit: fasta_pair
    tuple val(meta), path("*_interleaved.fasta.gz"), optional: true, emit: fasta_interleaved
    tuple val(meta), path("*_singleton.fasta.gz"), optional: true, emit: fasta_singleton
    tuple val(meta), path("*_other.fasta.gz"), optional: true, emit: fasta_other

    tuple val(meta), path("*_{1,2}.fastq.gz"), optional: true, emit: fastq_pair
    tuple val(meta), path("*_interleaved.fastq"), optional: true, emit: fastq_interleaved
    tuple val(meta), path("*_singleton.fastq.gz"), optional: true, emit: fastq_singleton
    tuple val(meta), path("*_other.fastq.gz"), optional: true, emit: fastq_other

    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_end = meta?.single_end ?: false

    // Subcommand configuration map
    //
    // Each entry describes how a samtools subcommand operates in terms of inputs and outputs.
    //
    // Keys:
    //
    //   stdout             Output flag used when the command is NOT the last in the command chain,
    //                      so that its output streams to the next command. '' if no flag is needed.
    //
    //   out                Flag placed before the output path when the command IS the last in
    //                      the chain. '' if no flag is needed (the path is positional, or unused).
    //
    //   threads            The command can accept the -@ threading argument.
    //
    //   uncompressed       Add "-u" (uncompressed output) when the command is not the last.
    //
    //   output_format_type Output format type, either 'alignment' (BAM/CRAM/SAM) or 'sequence' (FASTA/FASTQ).
    //
    def SUBCOMMAND_CONFIG = [
        view: [stdout: '', out: '-o', threads: true, uncompressed: true, output_format_type: 'alignment'],
        sort: [stdout: '', out: '-o', threads: true, uncompressed: true, output_format_type: 'alignment'],
        markdup: [stdout: '-', out: '', threads: true, uncompressed: true, output_format_type: 'alignment'],
        fixmate: [stdout: '-', out: '', threads: true, uncompressed: true, output_format_type: 'alignment'],
        merge: [stdout: '-o -', out: '-o', threads: true, uncompressed: true, output_format_type: 'alignment'],
        cat: [stdout: '', out: '-o', threads: false, uncompressed: false, output_format_type: 'alignment'],
        collate: [stdout: '-O', out: '-o', threads: true, uncompressed: true, output_format_type: 'alignment'],
        fasta: [stdout: '', out: '', threads: true, uncompressed: true, output_format_type: 'sequence'],
        fastq: [stdout: '', out: '', threads: true, uncompressed: true, output_format_type: 'sequence'],
    ]

    // Only allow chains with at least two commands
    def n_commands = chain.size()
    if (n_commands <= 1) {
        error("Error: SAMTOOLS_MULTICOMMAND requires at least two commands!")
    }

    // Validate chain commands
    // Check that all tools requested are allowed and that fasta/fastq fall at the end of the chain only
    chain
        .withIndex()
        .each { subcommand, idx ->
            if (!SUBCOMMAND_CONFIG.containsKey(subcommand)) {
                error("Error: ${subcommand} not a valid argument for SAMTOOLS_MULTICOMMAND! Valid options are: ${SUBCOMMAND_CONFIG.keySet().join(", ")}")
            }
            if (subcommand in ['merge', 'cat'] && idx != 0) {
                error("Error: ${subcommand} must be the first command for SAMTOOLS_MULTICOMMAND!")
            }
            if (SUBCOMMAND_CONFIG[subcommand].output_format_type == 'sequence' && idx != n_commands - 1) {
                error("Error: ${subcommand} must be the last command for SAMTOOLS_MULTICOMMAND!")
            }
        }

    // Extension arguments for each command: args for the first, args2 ... argsN for the rest
    def all_args = (0..<n_commands).collect { idx -> task.ext[idx == 0 ? "args" : "args${idx + 1}"] ?: "" }

    // Output format comes from --output-fmt in the final command's args, defaulting to BAM
    def output_extension = (all_args[-1] =~ /(?i)(?:--output-fmt[ =]|-O\s*)(sam|bam|cram)/).with { m -> m.find() ? m.group(1).toLowerCase() : "bam" }
    def is_cram_input = fasta && input.any { f -> f.getExtension() == "cram" }
    def is_cram_output = fasta && output_extension == "cram"

    def staged_names = [input].flatten().collect { f -> f.name }
    if (SUBCOMMAND_CONFIG[chain[-1]].output_format_type == 'alignment' && "${prefix}.${output_extension}" in staged_names) {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    // Output flag for a command, depending on whether it is last in the chain
    def output_option = { subcommand, cfg, is_last_command ->
        if (!is_last_command) {
            return cfg.stdout
        }
        if (cfg.output_format_type == 'sequence') {
            def seq_ext = "${subcommand}.gz"
            return [
                "-0 ${prefix}_other.${seq_ext}",
                "-1 ${prefix}_1.${seq_ext}",
                single_end ? "" : "-2 ${prefix}_2.${seq_ext}",
                "-s ${prefix}_singleton.${seq_ext}",
            ].findAll().join(" ")
        }
        return "${cfg.out} ${prefix}.${output_extension}".trim()
    }

    // Build the chained command string
    //
    // Options that differ per subcommand come from SUBCOMMAND_CONFIG. The rest depends on position
    // in the chain: the first command reads the input, the last writes the output, and
    // CRAM references are added where needed.
    //
    // EXAMPLE:
    //
    // This module will construct a chained samtools command from an input list of
    // subtools, such as [view, sort, markdup]:
    //
    //     samtools view ${args} input.bam |\
    //     samtools sort ${args2} |\
    //     samtools markdup ${args3} -o output.bam
    //
    // The args are numbered sequentially for each tool in the sequence and CRAM references
    // are automatically applied if needed. FASTA and FASTQ outputs are also available.
    def chained_command = chain
        .withIndex()
        .collect { subcommand, idx ->
            def cfg = SUBCOMMAND_CONFIG[subcommand]
            def cmd_args = all_args[idx]
            def is_first_command = idx == 0
            def is_last_command = idx == n_commands - 1

            def threads = (cfg.threads && !cmd_args.contains("-@")) ? "-@ ${task.cpus}" : ""
            def uncompressed = (!is_last_command && cfg.uncompressed) ? "-u" : ""
            def reference = ((is_first_command && is_cram_input) || (is_last_command && is_cram_output)) ? "--reference ${fasta}" : ""
            def cmd_input = is_first_command ? "${input}" : "-"
            def cmd_output = output_option.call(subcommand, cfg, is_last_command)

            return ["samtools", subcommand, threads, uncompressed, cmd_args, reference, cmd_input, cmd_output].findAll().join(" ")
        }
        .join(" |\\\n")

    """
    set -o pipefail
    ${chained_command}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single_end = meta?.single_end ?: false
    def n_commands = chain.size()
    def all_args = (0..<n_commands).collect { idx -> task.ext[idx == 0 ? "args" : "args${idx + 1}"] ?: "" }
    def output_extension = (all_args[-1] =~ /(?i)(?:--output-fmt[ =]|-O\s*)(sam|bam|cram)/).with { m -> m.find() ? m.group(1).toLowerCase() : "bam" }
    def final_command = chain[-1]

    def stub_outputs = []

    if (final_command in ['fasta', 'fastq']) {
        def ext = final_command == 'fasta' ? 'fasta.gz' : 'fastq.gz'
        stub_outputs << "echo | bgzip > ${prefix}_1.${ext}"
        if (!single_end) {
            stub_outputs << "echo | bgzip > ${prefix}_2.${ext}"
        }
        stub_outputs << "echo | bgzip > ${prefix}_singleton.${ext}"
        stub_outputs << "echo | bgzip > ${prefix}_other.${ext}"
    }
    else {
        stub_outputs << "touch ${prefix}.${output_extension}"
    }
    """
    ${stub_outputs.join("\n")}
    """
}
