// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MAXBIN2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::maxbin2=2.2.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/maxbin2:2.2.7--he1b5a44_2"
    } else {
        container "quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_2"
    }

    input:
    tuple val(meta), path(contigs), path(reads), path(abund)

    output:
    tuple val(meta), path("*.fasta.gz")   , emit: binned_fastas
    tuple val(meta), path("*.summary")    , emit: summary
    tuple val(meta), path("*.log.gz")     , emit: log
    tuple val(meta), path("*.marker.gz")  , emit: marker_counts
    tuple val(meta), path("*.noclass.gz") , emit: unbinned_fasta
    tuple val(meta), path("*.tooshort.gz"), emit: tooshort_fasta
    tuple val(meta), path("*_bin.tar.gz") , emit: marker_bins , optional: true
    tuple val(meta), path("*_gene.tar.gz"), emit: marker_genes, optional: true
    path "versions.yml"                   , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def associate_files = reads ? "-reads $reads" : "-abund $abund"
    """
    run_MaxBin.pl \\
        -contig $contigs \\
        $associate_files \\
        -thread $task.cpus \\
        $options.args \\
        -out $prefix

    gzip *.fasta *.noclass *.tooshort *log *.marker

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        maxbin2: \$( run_MaxBin.pl -v | head -n 1 | sed 's/MaxBin //' )
    END_VERSIONS
    """
}
