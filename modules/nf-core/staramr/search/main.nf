process STARAMR_SEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::staramr=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staramr:0.9.1--pyhdfd78af_0':
        'biocontainers/staramr:0.9.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly_fasta) // assembly as fasta file

    output:
    tuple val(meta), path("out/results.xlsx")        , emit: xlsx
    tuple val(meta), path("out/summary.tsv")         , emit: summary_tsv
    tuple val(meta), path("out/detailed_summary.tsv"), emit: detailed_summary_tsv
    tuple val(meta), path("out/resfinder.tsv")       , emit: resfinder_tsv
    tuple val(meta), path("out/plasmidfinder.tsv")   , emit: plasmidfinder_tsv
    tuple val(meta), path("out/mlst.tsv")            , emit: mlst_tsv
    tuple val(meta), path("out/settings.txt")        , emit: txt
    tuple val(meta), path("out/pointfinder.tsv")     , emit: pointfinder_tsv         , optional: true
    tuple val(meta), path("out/hits/resfinder*")     , emit: hits_resfinder_fasta    , optional: true
    tuple val(meta), path("out/hits/pointfinder*")   , emit: hits_pointfinder_fasta  , optional: true
    tuple val(meta), path("out/hits/plasmidfinder*") , emit: hits_plasmidfinder_fasta, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${assembly_fasta} > ${assembly_fasta}.uncompressed
    staramr \\
        search \\
        $args \\
        --nprocs $task.cpus \\
        -o out \\
        ${assembly_fasta}.uncompressed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staramr : \$(echo \$(staramr --version 2>&1) | sed 's/^.*staramr //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir out
    touch out/results.xlsx
    touch out/{summary,detailed_summary,resfinder,pointfinder,plasmidfinder,mlst}.tsv
    touch settings.txt
    mkdir out/hits

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staramr : \$(echo \$(staramr --version 2>&1) | sed 's/^.*staramr //' )
    END_VERSIONS
    """
}
