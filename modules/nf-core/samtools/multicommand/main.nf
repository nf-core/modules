process SAMTOOLS_MULTICOMMAND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(input, arity: '1..*'), path(index, arity: '0..*')
    tuple val(meta2), path(fasta), path(fai)
    val pipeline

    output:
    // Alignment format outputs (view, sort, markdup, merge, cat, collate)
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.cram"), optional: true, emit: cram
    tuple val(meta), path("*.sam"), optional: true, emit: sam
    tuple val(meta), path("*.{bai,csi,crai}"), optional: true, emit: index

    // Sequence outputs (fasta, fastq)
    tuple val(meta), path("*.fasta.gz"), emit: fasta, optional: true
    tuple val(meta), path("*.fastq.gz"), emit: fastq, optional: true
    tuple val(meta), path("*_{1,2}.fasta.gz"), emit: fasta_pair, optional: true
    tuple val(meta), path("*_interleaved.fasta.gz"), emit: fasta_interleaved, optional: true
    tuple val(meta), path("*_singleton.fasta.gz"), emit: fasta_singleton, optional: true
    tuple val(meta), path("*_other.fasta.gz"), emit: fasta_other, optional: true
    tuple val(meta), path("*_{1,2}.fastq.gz"), emit: fastq_pair, optional: true
    tuple val(meta), path("*_interleaved.fastq"), emit: fastq_interleaved, optional: true
    tuple val(meta), path("*_singleton.fastq.gz"), emit: fastq_singleton, optional: true
    tuple val(meta), path("*_other.fastq.gz"), emit: fastq_other, optional: true

    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (pipeline.size() <= 1) {
        error("Error: SAMTOOLS_MULTICOMMAND requires at least two commands!")
    }

    def valid_options = ['view', 'sort', 'markdup', 'fixmate', 'merge', 'cat', 'collate', 'fastq', 'fasta']
    pipeline.collect { tool ->
        if (!(tool in valid_options)) {
            error("Error: ${tool} not a valid pipeline argument for SAMTOOLS_MULTICOMMAND! Valid options are: ${valid_options.join(", ")}")
        }
    }

    def n_commands = pipeline.size()
    def is_cram_input = fasta && input.collect { f -> f.getExtension() == "cram" }.any()
    def is_cram_output = fasta && get_args(task, n_commands - 1)

    // Build the pipeline command
    //
    // Note: The functions used to compile the pipeline are defined below the process definition.
    // They have been separated into different concerns in order to effectively capture exceptions
    // to the normal operation of samtools commands - for example, samtools collate uses
    // "-O" to write to stdout rather than "-". This makes capturing further exceptions easier to manage
    // without requiring complex conditional logic.
    //
    def pipeline_command = pipeline
        .withIndex()
        .collect { subcommand, idx ->
            def is_first_command = (idx == 0)
            def is_last_command = (idx == n_commands - 1)

            def cmd_args = get_args(task.ext, idx)
            def cmd_threads = get_threads(subcommand, task.cpus, args)
            def cmd_input = is_first_command ? "${input}" : get_stdin()
            def cmd_output = is_last_command ? get_file_output(subcommand, cmd_args, prefix, meta?.single_end ?: false) : get_stdout(subcommand)
            def uncompressed = !is_last_command ? get_uncompressed_flag(subcommand) : ""
            def reference = get_reference(is_first_command, is_last_command, is_cram_input, is_cram_output, fasta)

            return ["samtools", subcommand, cmd_threads, uncompressed, cmd_args, reference, cmd_input, cmd_output].findAll().join(" ")
        }
        .join(" |\\\n")

    // EXAMPLE:
    //
    // This module will construct a samtools pipeline command from an input list of
    // subtools, such as [view, sort, markdup]:
    //
    //     samtools view ${args} input.bam |\
    //     samtools sort ${args2} |\
    //     samtools markdup ${args3} -o output.bam
    //
    // The args are numbered sequentially for each tool in the sequence and CRAM references
    // are automatically applied if needed. FASTA and FASTQ outputs are also available.
    """
    ${pipeline_command}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def valid_options = ['view', 'sort', 'markdup', 'fixmate', 'merge', 'cat', 'collate', 'fastq', 'fasta']
    pipeline.collect { tool ->
        if (!(tool in valid_options)) {
            error("Error: ${tool} not a valid pipeline argument for SAMTOOLS_PIPELINE! Valid options are: ${valid_options.join(", ")}")
        }
    }

    def n_commands = pipeline.size()
    def final_command = pipeline[n_commands - 1]
    def final_args = get_args(task, n_commands - 1)

    def stub_outputs = []

    if (final_command in ['view', 'sort', 'merge', 'cat', 'markdup', 'fixmate', 'collate']) {
        def extension = get_output_extension(final_args)
        stub_outputs << "touch ${prefix}.${extension}"
    }
    else if (final_command == "fasta") {
        if (meta.single_end) {
            stub_outputs << "echo | gzip > ${prefix}_1.fasta.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fasta.gz"
        }
        else {
            stub_outputs << "echo | gzip > ${prefix}_1.fasta.gz"
            stub_outputs << "echo | gzip > ${prefix}_2.fasta.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fasta.gz"
        }
        stub_outputs << "echo | gzip > ${prefix}_other.fasta.gz"
    }
    else if (final_command == "fastq") {
        if (meta.single_end) {
            stub_outputs << "echo | gzip > ${prefix}_1.fastq.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fastq.gz"
        }
        else {
            stub_outputs << "echo | gzip > ${prefix}_1.fastq.gz"
            stub_outputs << "echo | gzip > ${prefix}_2.fastq.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fastq.gz"
        }
        stub_outputs << "echo | gzip > ${prefix}_other.fastq.gz"
    }

    """
    ${stub_outputs.join("\n")}
    """
}

// ============================================================================
// Helper Functions
// ============================================================================


// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Returns stdin redirection symbol for piping between samtools commands.
 * @return String "-" for stdin
 */
def get_stdin() {
    return "-"
}

/**
 * Determines stdout redirection flag for a samtools subcommand.
 * @param subcommand The samtools subcommand (e.g. 'view', 'sort', 'collate')
 * @return String Redirection flag: "-O" for collate, "" for view/sort/cat, "-" otherwise
 */
def get_stdout(subcommand) {
    if (subcommand == "collate") {
        return "-O"
    }
    else if (subcommand in ["view", "sort", "cat"]) {
        return ""
    }
    return "-"
}

/**
 * Determines output file format based on task arguments.
 * @param args Task extension arguments string
 * @return String Output format: 'sam', 'cram', or 'bam' (default)
 */
def get_output_extension(args) {
    if (args.contains("--output-fmt sam")) {
        return "sam"
    }
    else if (args.contains("--output-fmt cram")) {
        return "cram"
    }
    return "bam"
}

/**
 * Generates output file arguments for a samtools subcommand.
 * @param subcommand The samtools subcommand
 * @param args Task extension arguments string
 * @param prefix Output file prefix
 * @param single_end Whether input is single-end
 * @return String Output file argument(s) for the subcommand
 */
def get_file_output(subcommand, args, prefix, single_end) {
    def extension = get_output_extension(args)

    if (subcommand in ['view', 'sort', 'merge', 'cat', 'fixmate', 'collate']) {
        return "-o ${prefix}.${extension}"
    }
    else if (subcommand == "markdup") {
        return "${prefix}.${extension}"
    }
    else if (subcommand == "fasta") {
        def output = "-0 ${prefix}_other.fasta.gz"
        if (!single_end) {
            output += " -1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz -s ${prefix}_singleton.fasta.gz"
        }
        else {
            output += " -1 ${prefix}_1.fasta.gz -s ${prefix}_singleton.fasta.gz"
        }
        return output
    }
    else if (subcommand == "fastq") {
        def output = "-0 ${prefix}_other.fastq.gz"
        if (!single_end) {
            output += " -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
        }
        else {
            output += " -1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz"
        }
        return output
    }

    error("Trying to get output for unsupported subcommand: ${subcommand}")
}

/**
 * Retrieves task extension arguments for a given pipeline step index.
 * @param task Nextflow task object
 * @param idx Zero-based pipeline step index
 * @return String Task extension arguments, or empty string if not defined
 */
def get_args(ext, idx) {
    def argsKey = idx == 0 ? "args" : "args${idx + 1}"
    return ext[argsKey] ?: ""
}

/**
 * Generates thread argument for a samtools subcommand.
 * @param subcommand The samtools subcommand
 * @param task Nextflow task object
 * @return String Thread argument (e.g. "-@ 4"), or empty string for cat
 */
def get_threads(subcommand, cpus, args) {
    if (subcommand != "cat") {
        if (!args.contains("-@")) {
            return "-@ ${cpus}"
        }
    }
    return ""
}

/**
 * Determines whether to add uncompressed output flag for intermediate commands.
 * @param subcommand The samtools subcommand
 * @return String Uncompressed flag "-u", or empty string for cat
 */
def get_uncompressed_flag(subcommand) {
    if (subcommand != "cat") {
        return "-u"
    }
    return ""
}

/**
 * Generates reference FASTA argument if needed for CRAM input or output.
 * @param is_first_command Whether this is the first command in the pipeline
 * @param is_last_command Whether this is the last command in the pipeline
 * @param is_cram_input Whether the input is CRAM format
 * @param is_cram_output Whether the output is CRAM format
 * @param fasta Path to reference FASTA file
 * @return String Reference argument (e.g. "--reference file.fa"), or empty string if not needed
 */
def get_reference(is_first_command, is_last_command, is_cram_input, is_cram_output, fasta) {
    def reference_string = "--reference ${fasta}"
    if (is_cram_input && is_first_command) {
        return reference_string
    }
    else if (is_cram_output && is_last_command) {
        return reference_string
    }
    return ""
}
