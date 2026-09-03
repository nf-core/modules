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
    def final_command = pipeline[n_commands - 1]
    def is_cram_input = fasta && input.collect { f -> f.getExtension() == "cram" }.any()

    // Build the pipeline command
    def pipeline_command = pipeline
        .withIndex()
        .collect { subcommand, idx ->
            def is_first = (idx == 0)
            def is_last = (idx == n_commands - 1)
            def args = get_args(task, idx)
            def input_part = get_input(is_first, input)
            def output_part = get_output(subcommand, is_last, final_command, task, n_commands, prefix, meta?.single_end ?: false)
            def threads_part = get_threads(subcommand, task)
            def reference_part = get_reference(is_first, is_last, is_cram_input, args, fasta)
            def uncompressed = !is_last && subcommand != "cat" ? "-u" : ""

            return build_command(subcommand, threads_part, args, input_part, output_part, reference_part, uncompressed)
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

    def stub_outputs = []

    if (final_command in ['view', 'sort', 'merge', 'cat', 'markdup', 'fixmate', 'collate']) {
        def argsKey = n_commands == 1 ? "args" : "args${n_commands}"
        def argsLast = task.ext[argsKey] ?: ""
        def extension = argsLast.contains("--output-fmt sam")
            ? "sam"
            : argsLast.contains("--output-fmt cram")
                ? "cram"
                : "bam"
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

def get_input(is_first, input) {
    if (is_first) {
        return input.join(" ")
    }
    return "-"
}

def get_output(subcommand, is_last, final_command, task, n_commands, prefix, single_end) {
    if (!is_last) {
        if (subcommand == "collate") {
            return "-O"
        }
        else if (subcommand in ["view", "sort", "cat"]) {
            return ""
        }
        return "-"
    }

    def argsKey = n_commands == 1 ? "args" : "args${n_commands}"
    def argsLast = task.ext[argsKey] ?: ""

    if (final_command in ['view', 'sort', 'merge', 'cat', 'fixmate', 'collate']) {
        def extension = argsLast.contains("--output-fmt sam")
            ? "sam"
            : argsLast.contains("--output-fmt cram")
                ? "cram"
                : "bam"
        return "-o ${prefix}.${extension}"
    }
    else if (final_command == "markdup") {
        def extension = argsLast.contains("--output-fmt sam")
            ? "sam"
            : argsLast.contains("--output-fmt cram")
                ? "cram"
                : "bam"
        return "${prefix}.${extension}"
    }
    else if (final_command == "fasta") {
        def output = "-0 ${prefix}_other.fasta.gz"
        if (!single_end) {
            output += " -1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz -s ${prefix}_singleton.fasta.gz"
        }
        else {
            output += " -1 ${prefix}_1.fasta.gz -s ${prefix}_singleton.fasta.gz"
        }
        return output
    }
    else if (final_command == "fastq") {
        def output = "-0 ${prefix}_other.fastq.gz"
        if (!single_end) {
            output += " -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
        }
        else {
            output += " -1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz"
        }
        return output
    }

    return ""
}

def get_args(task, idx) {
    def argsKey = idx == 0 ? "args" : "args${idx + 1}"
    return task.ext[argsKey] ?: ""
}

def get_threads(subcommand, task) {
    return subcommand != "cat" ? "-@ ${task.cpus}" : ""
}

def get_reference(is_first, is_last, is_cram_input, args, fasta) {
    if (is_first && is_cram_input) {
        return "--reference ${fasta}"
    }
    if (is_last && args.contains("--output-fmt cram")) {
        return "--reference ${fasta}"
    }
    return ""
}

def build_command(subcommand, threads, args, input, output, reference, uncompressed) {
    def parts = ["samtools", subcommand]

    if (threads) {
        parts << threads
    }
    if (args) {
        parts << args
    }
    if (input) {
        parts << input
    }
    if (reference) {
        parts << reference
    }
    if (output) {
        parts << output
    }
    if (uncompressed) {
        parts << uncompressed
    }

    return parts.findAll().join(" ")
}
