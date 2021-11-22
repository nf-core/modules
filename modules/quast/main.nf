process QUAST {
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::quast=5.0.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2'
    } else {
        container 'quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2'
    }

    input:
    path consensus
    path fasta
    path gff
    val use_fasta
    val use_gff

    output:
    path "${prefix}"    , emit: results
    path '*.tsv'        , emit: tsv
    path "versions.yml" , emit: versions

    script:
    prefix        = options.suffix ?: software
    def features  = use_gff ? "--features $gff" : ''
    def reference = use_fasta ? "-r $fasta" : ''
    """
    quast.py \\
        --output-dir $prefix \\
        $reference \\
        $features \\
        --threads $task.cpus \\
        $args \\
        ${consensus.join(' ')}
    ln -s ${prefix}/report.tsv
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
