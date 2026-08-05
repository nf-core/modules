process PARSESDRF_CONVERT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sdrf-pipelines:0.1.4--pyhdfd78af_0':
        'quay.io/biocontainers/sdrf-pipelines:0.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sdrf), path(fasta), val(raw_folder)
    val format

    output:
    tuple val(meta), path("${prefix}_samplesheet.tsv"), path("${prefix}_presets.tsv"),             emit: mhcquant,     optional: true
    tuple val(meta), path("${prefix}_samplesheet.tsv"), path("${prefix}_experimental_design.tsv"), emit: openms,       optional: true
    tuple val(meta), path("${prefix}.xml"),             path("${prefix}_design.txt"),              emit: maxquant,     optional: true
    tuple val(meta), path("${prefix}.csv"),                                                        emit: msstats,      optional: true
    tuple val(meta), path("${prefix}_design.csv"),      path("${prefix}_comparisons.csv"),         emit: normalyzerde, optional: true
    tuple val(meta), path("${prefix}.cfg"),             path("${prefix}_design.tsv"),              emit: diann,        optional: true
    tuple val("${task.process}"), val('sdrf-pipelines'), eval("parse_sdrf --version | cut -d ' ' -f 2"), topic: versions, emit: versions_sdrfpipelines

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def output_args = ''
    def post = ''

    if (format == 'mhcquant') {
        output_args = "-os ${prefix}_samplesheet.tsv -op ${prefix}_presets.tsv"
    } else if (format == 'openms') {
        // convert-openms writes openms.tsv + (optionally) experimental_design.tsv
        post = "mv openms.tsv ${prefix}_samplesheet.tsv ; [ -f experimental_design.tsv ] && mv experimental_design.tsv ${prefix}_experimental_design.tsv || true"
    } else if (format == 'maxquant') {
        if (!fasta || !raw_folder) {
            error "PARSESDRF_CONVERT: format 'maxquant' requires both `fasta` (path) and `raw_folder` (string) inputs"
        }
        // fasta is staged into the work dir; raw_folder is a literal string embedded into mqpar.xml's <filePaths>
        output_args = "-f ${fasta} -r ${raw_folder} -o1 ${prefix}.xml -o2 ${prefix}_design.txt"
    } else if (format == 'msstats') {
        output_args = "-o ${prefix}.csv"
    } else if (format == 'normalyzerde') {
        output_args = "-o ${prefix}_design.csv -oc ${prefix}_comparisons.csv"
    } else if (format == 'diann') {
        // convert-diann writes diann_config.cfg + diann_design.tsv with fixed names
        post = "mv diann_config.cfg ${prefix}.cfg && mv diann_design.tsv ${prefix}_design.tsv"
    } else {
        error "PARSESDRF_CONVERT: unsupported format '${format}' (expected one of: mhcquant, openms, maxquant, msstats, normalyzerde, diann)"
    }

    """
    parse_sdrf convert-${format} \\
        -s ${sdrf} \\
        ${output_args} \\
        ${args}
    ${post}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def stub_files
    if (format == 'mhcquant') {
        stub_files = ["${prefix}_samplesheet.tsv", "${prefix}_presets.tsv"]
    } else if (format == 'openms') {
        stub_files = ["${prefix}_samplesheet.tsv", "${prefix}_experimental_design.tsv"]
    } else if (format == 'maxquant') {
        stub_files = ["${prefix}.xml", "${prefix}_design.txt"]
    } else if (format == 'msstats') {
        stub_files = ["${prefix}.csv"]
    } else if (format == 'normalyzerde') {
        stub_files = ["${prefix}_design.csv", "${prefix}_comparisons.csv"]
    } else if (format == 'diann') {
        stub_files = ["${prefix}.cfg", "${prefix}_design.tsv"]
    } else {
        error "PARSESDRF_CONVERT: unsupported format '${format}'"
    }
    def touches = stub_files.collect { f -> "touch ${f}" }.join('\n    ')
    """
    ${touches}
    """
}
