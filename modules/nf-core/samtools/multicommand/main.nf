process SAMTOOLS_MULTICOMMAND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta), path(fai)
    val(pipeline)

    output:
    // Alignment format outputs (view, sort, markdup, merge, cat, collate)
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.cram"), optional: true, emit: cram
    tuple val(meta), path("*.sam"), optional: true, emit: sam
    tuple val(meta), path("*.{bai,csi,crai}"), optional: true, emit: index

    // Sequence outputs (fasta, fastq)
    tuple val(meta), path("*.fasta.gz"), optional: true, emit: fasta
    tuple val(meta), path("*.fastq.gz"), optional: true, emit: fastq

    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def valid_options = ['view', 'sort', 'markdup', 'fixmate', 'merge', 'cat', 'collate', 'fastq', 'fasta']
    pipeline.collect { tool ->
        if (!(tool in valid_options)) {
            error("Error: ${tool} not a valid pipeline argument for SAMTOOLS_PIPELINE! Valid options are: ${valid_options.join(", ")}")
        }
    }

    def n_commands = pipeline.size()
    def final_command = pipeline[n_commands - 1]

    // Build output string based on final command
    def output_string = ""
    def input_reference = (fasta && input.getExtension() == "cram") ? "--reference ${fasta}" : ""
    def output_reference = ""

    if (final_command in ['view', 'sort', 'merge', 'cat', 'markdup', 'fixmate', 'merge', 'cat', 'collate']) {
        // These produce alignment files
        def argsKey = n_commands == 1 ? "args" : "args${n_commands}"
        def argsLast = task.ext[argsKey] ?: ""
        def extension = argsLast.contains("--output-fmt sam")
            ? "sam"
            : argsLast.contains("--output-fmt cram")
                ? "cram"
                : "bam"
        output_reference = (fasta && input.getExtension() == "cram") ? "--reference ${fasta}" : ""
        output_string = "-o ${prefix}.${extension}"
    } else if (final_command == "fasta") {
        // fasta produces multiple files with special output flags
        output_string = "-0 ${prefix}_other.fasta.gz"
        if (!meta.single_end) {
            output_string = output_string + " -1 ${prefix}_1.fasta.gz -2 ${prefix}_2.fasta.gz -s ${prefix}_singleton.fasta.gz"
        } else {
            output_string = output_string + " -1 ${prefix}_1.fasta.gz -s ${prefix}_singleton.fasta.gz"
        }
    } else if (final_command == "fastq") {
        // fastq produces multiple files with special output flags
        output_string = "-0 ${prefix}_other.fastq.gz"
        if (!meta.single_end) {
            output_string = output_string + " -1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
        } else {
            output_string = output_string + " -1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz"
        }
    }

    // Build the pipeline command
    def pipeline_command = pipeline.withIndex().collect { subcommand, idx ->
        def argsKey = idx == 0 ? "args" : "args${idx + 1}"
        def taskArgs = task.ext[argsKey] ?: ""

        def cmd_parts = ["samtools", subcommand]
        if (taskArgs) {
            cmd_parts << taskArgs
        }
        if (idx == 0) {
            if (input_reference) {
                cmd_parts << input_reference
            }
            cmd_parts << (input instanceof List ? input.join(" ") : input)
        }
        if (idx == n_commands - 1) {
            if (output_reference) {
                cmd_parts << output_reference
            }
            cmd_parts << output_string
        }

        return cmd_parts.join(" ")
    }.join(" |\\\n")

    // EXAMPLE:
    //
    // This module will construct a samtools pipeline command from an input list of
    // subtools, such as [view, sort, markdup]:
    //
    //     samtools view ${args} input.bam |\
    //     samtools sort ${args2} |\
    //     samtools markdup ${args3} -o output.bam
    //
    // The args are numbered sequenctially for each tool in the sequence and CRAM references
    // are automatically applied if needed. FASTA and FASTQ outputs are also available.
    """
    ${pipeline_command}
    """

    stub:
    def args = task.ext.args ?: ''
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
    } else if (final_command == "fasta") {
        if (meta.single_end) {
            stub_outputs << "echo | gzip > ${prefix}_1.fasta.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fasta.gz"
        } else {
            stub_outputs << "echo | gzip > ${prefix}_1.fasta.gz"
            stub_outputs << "echo | gzip > ${prefix}_2.fasta.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fasta.gz"
        }
        stub_outputs << "echo | gzip > ${prefix}_other.fasta.gz"
    } else if (final_command == "fastq") {
        if (meta.single_end) {
            stub_outputs << "echo | gzip > ${prefix}_1.fastq.gz"
            stub_outputs << "echo | gzip > ${prefix}_singleton.fastq.gz"
        } else {
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
