// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// TOOL DESCRIPTION:
// Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It was developed for the use in a clinical research setting.
// Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, such as viral integration sites,
// internal tandem duplications, whole exon duplications, truncations of genes (i.e., breakpoints in introns and intergenic regions).
// Source: https://github.com/suhrig/arriba
// Article: Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research.
// DOI: 10.1101/gr.257246.119

params.options = [:]
options        = initOptions(params.options)

process ARRIBA {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::arriba=2.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/arriba:2.1.0--h3198e80_1"
    } else {
        container "quay.io/biocontainers/arriba:2.1.0--h3198e80_1"
    }

    input:
    tuple val(meta), path(bam)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*_fusions.tsv")       , emit: tsv
    path "*.version.txt"                         , emit: version

    script:
    def software            = getSoftwareName(task.process)
    def prefix              = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def blacklist_option    = (options.args.contains('-b')) ? '' : '-f blacklist'

    """
    arriba \\
        -x $bam \\
        -a $fasta \\
        -g $gtf \\
        -o ${prefix}_fusions.tsv \\
        $blacklist_option \\
        $options.args

    echo \$(arriba -h | grep 'Version:' 2>&1) |  sed 's/Version:\s//' > ${software}.version.txt
    """
}
