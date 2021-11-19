include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DASTOOL_DASTOOL {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::das_tool=1.1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/das_tool:1.1.3--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/das_tool:1.1.3--r41hdfd78af_0"
    }

    input:
    tuple val(meta), path(contigs), path(bins)
    path(proteins)
    path(db_directory)
    val(search_engine)

    output:
    tuple val(meta), path("*.log")                                      , emit: log
    tuple val(meta), path("*_summary.txt")                              , emit: summary
    tuple val(meta), path("*_DASTool_scaffolds2bin.txt")                , emit: scaffolds2bin
    tuple val(meta), path("*.eval")                     , optional: true, emit: eval
    tuple val(meta), path("*_DASTool_bins/*.fa")        , optional: true, emit: bins
    tuple val(meta), path("*.pdf")                      , optional: true, emit: pdfs
    tuple val(meta), path("*.proteins.faa")             , optional: true, emit: fasta_proteins
    tuple val(meta), path("*.archaea.scg")              , optional: true, emit: fasta_archaea_scg
    tuple val(meta), path("*.bacteria.scg")             , optional: true, emit: fasta_bacteria_scg
    tuple val(meta), path("*.seqlength")                , optional: true, emit: seqlength
    path "versions.yml"                                 , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bin_list = bins instanceof List ? bins.join(",") : "$bins"
    def engine = search_engine ? "--search_engine $search_engine" : "--search_engine diamond"
    def db_dir = db_directory ? "--db_directory $db_directory" : ""
    def clean_contigs = contigs.toString() - ".gz"
    def decompress_contigs = contigs.toString() == clean_contigs ? "" : "gunzip -q -f $contigs"
    def decompress_proteins = proteins ? "gunzip -f $proteins" : ""
    def clean_proteins = proteins ? proteins.toString() - ".gz" : ""
    def proteins_pred = proteins ? "--proteins $clean_proteins" : ""

    if (! search_engine) {
        log.info('[DAS_Tool] Default search engine (USEARCH) is proprietary software and not available in bioconda. Using DIAMOND as alternative.')
    }

    """
    $decompress_proteins
    $decompress_contigs

    DAS_Tool \\
        $options.args \\
        $proteins_pred \\
        $db_dir \\
        $engine \\
        -t $task.cpus \\
        --bins $bin_list \\
        -c $clean_contigs \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( DAS_Tool --version 2>&1 | grep "DAS Tool" | sed 's/DAS Tool version //' )
    END_VERSIONS
    """
}
