process SAMTOOLS_PIPELINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam)
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
    tuple val(meta), path("*_interleaved.*"), optional: true, emit: interleaved

    // Statistics outputs (coverage, depth, stats, idxstats, flagstat)
    tuple val(meta), path("*.txt"), optional: true, emit: coverage
    tuple val(meta), path("*.tsv"), optional: true, emit: depth
    tuple val(meta), path("*.stats"), optional: true, emit: stats
    tuple val(meta), path("*.idxstats"), optional: true, emit: idxstats
    tuple val(meta), path("*.flagstat"), optional: true, emit: flagstat

    tuple val("${task.process}"), val('samtools'), eval("samtools --version"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_commands = pipeline.size()
    def final_command = pipeline[n_commands - 1]

    // Build output string based on final command
    def output_string = ""

    if (final_command in ['view', 'sort', 'markdup', 'fixmate', 'merge', 'cat', 'collate']) {
        // These produce alignment files
        def argsLast = task.ext["args${n_commands - 1}"] ?: ""
        def extension = argsLast.contains("--output-fmt sam")
            ? "sam"
            : argsLast.contains("--output-fmt cram")
                ? "cram"
                : "bam"
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
    } else if (final_command == "coverage") {
        // coverage produces a txt file
        output_string = "-o ${prefix}.txt"
    } else if (final_command == "depth") {
        // depth produces a tsv file
        output_string = "-o ${prefix}.tsv"
    } else if (final_command == "stats") {
        // stats produces a stats file
        output_string = "> ${prefix}.stats"
    } else if (final_command == "idxstats") {
        // idxstats produces an idxstats file
        output_string = "> ${prefix}.idxstats"
    } else if (final_command == "flagstat") {
        // flagstat produces a flagstat file
        output_string = "> ${prefix}.flagstat"
    }

    // Build the pipeline command
    def pipeline_command = pipeline.withIndex().collect { subcommand, idx ->
        def argsKey = "args${idx}"
        def taskArgs = task.ext[argsKey] ?: ""

        def cmd = ""
        if (idx == 0) {
            // First command: read from input BAM
            cmd = "samtools ${subcommand} ${taskArgs} ${bam}"
        } else {
            // Subsequent commands: read from stdin
            cmd = "samtools ${subcommand} ${taskArgs}"
        }

        // Last command: add output string
        if (idx == n_commands - 1) {
            cmd = cmd + " ${output_string}"
        }

        cmd
    }.join(" |\\\n")
    """
    ${pipeline_command}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_commands = pipeline.size()
    def final_command = pipeline[n_commands - 1]

    def stub_outputs = ""

    if (final_command in ['view', 'sort', 'markdup', 'fixmate', 'collate']) {
        def argsLast = task.ext["args${n_commands - 1}"] ?: ""
        def extension = argsLast.contains("--output-fmt sam")
            ? "sam"
            : argsLast.contains("--output-fmt cram")
                ? "cram"
                : "bam"
        stub_outputs = "touch ${prefix}.${extension}"
    } else if (final_command == "merge") {
        def file_type = bam instanceof List ? bam[0].getExtension() : bam.getExtension()
        stub_outputs = "touch ${prefix}.${file_type}"
    } else if (final_command == "cat") {
        def file_type = bam instanceof List ? bam[0].getExtension() : bam.getExtension()
        stub_outputs = "touch ${prefix}.${file_type}"
    } else if (final_command == "fasta") {
        if (meta.single_end) {
            stub_outputs = "echo > ${prefix}_1.fasta.gz"
            stub_outputs = "echo > ${prefix}_singleton.fasta.gz"
        } else {
            stub_outputs = "echo > ${prefix}_1.fasta.gz"
            stub_outputs = "echo > ${prefix}_2.fasta.gz"
            stub_outputs = "echo > ${prefix}_singleton.fasta.gz"
        }
    } else if (final_command == "fastq") {
        if (meta.single_end) {
            stub_outputs = "echo > ${prefix}_1.fastq.gz"
            stub_outputs = "echo > ${prefix}_singleton.fastq.gz"
        } else {
            stub_outputs = "echo > ${prefix}_1.fastq.gz"
            stub_outputs = "echo > ${prefix}_2.fastq.gz"
            stub_outputs = "echo > ${prefix}_singleton.fastq.gz"
        }
        stub_outputs = "touch ${prefix}_other.fastq.gz"
    } else if (final_command == "coverage") {
        stub_outputs = "touch ${prefix}.txt"
    } else if (final_command == "depth") {
        stub_outputs = "touch ${prefix}.tsv"
    } else if (final_command == "stats") {
        stub_outputs = "touch ${prefix}.stats"
    } else if (final_command == "idxstats") {
        stub_outputs = "touch ${prefix}.idxstats"
    } else if (final_command == "flagstat") {
        stub_outputs = "touch ${prefix}.flagstat"
    }

    """
    ${stub_outputs}
    """
}
