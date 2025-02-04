process MAXBIN2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/maxbin2:2.2.7--he1b5a44_2' :
        'biocontainers/maxbin2:2.2.7--he1b5a44_2' }"

    input:
    tuple val(meta), path(contigs), path(reads), path(abund)

    output:
    tuple val(meta), path("*.fasta.gz")   , emit: binned_fastas
    tuple val(meta), path("*.summary")    , emit: summary
    tuple val(meta), path("*.abundance")  , emit: abundance   , optional: true
    tuple val(meta), path("*.log.gz")     , emit: log
    tuple val(meta), path("*.marker.gz")  , emit: marker_counts
    tuple val(meta), path("*.noclass.gz") , emit: unbinned_fasta
    tuple val(meta), path("*.tooshort.gz"), emit: tooshort_fasta
    tuple val(meta), path("*_bin.tar.gz") , emit: marker_bins , optional: true
    tuple val(meta), path("*_gene.tar.gz"), emit: marker_genes, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (reads && abund) { error("ERROR: MaxBin2 can only accept one of `reads` or `abund`, no both. Check input.") }
    def associate_files = ""
    if (reads)  {
        associate_files = "-reads $reads"
    } else if (abund instanceof List) {
        associate_files = [0..(abund.size() - 1)].collect { n ->
            def arg_n = n == 0 ? "" : "${n}"
            return " -${abund}${arg_n} ${abund[n]}"
        }.join()
    } else {
        associate_files = "-abund $abund"
    }
    """
    mkdir input/ && mv $contigs input/
    run_MaxBin.pl \\
        -contig input/$contigs \\
        $associate_files \\
        -thread $task.cpus \\
        $args \\
        -out $prefix

    gzip *.fasta *.noclass *.tooshort *log *.marker

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maxbin2: \$( run_MaxBin.pl -v | head -n 1 | sed 's/MaxBin //' )
    END_VERSIONS
    """
}
