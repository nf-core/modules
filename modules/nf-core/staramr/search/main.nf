// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process STARAMR_SEARCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02718d049d5a1f303ec4ef1f889aa3c96b33924133150862bd65420da000eb2a/data'
        : 'community.wave.seqera.io/library/mlst_staramr_biopython_numpy_pandas:f28fd6fb415e425b'}"

    input:
    tuple val(meta), path(genome_fasta) // genome as a fasta file

    output:
    tuple val(meta), path("*_results/results.xlsx")        , emit: results_xlsx
    tuple val(meta), path("*_results/summary.tsv")         , emit: summary_tsv
    tuple val(meta), path("*_results/detailed_summary.tsv"), emit: detailed_summary_tsv
    tuple val(meta), path("*_results/resfinder.tsv")       , emit: resfinder_tsv
    tuple val(meta), path("*_results/plasmidfinder.tsv")   , emit: plasmidfinder_tsv
    tuple val(meta), path("*_results/mlst.tsv")            , emit: mlst_tsv
    tuple val(meta), path("*_results/settings.txt")        , emit: settings_txt
    tuple val(meta), path("*_results/pointfinder.tsv")     , emit: pointfinder_tsv, optional: true
    tuple val("${task.process}"), val('staramr'), eval("staramr --version 2>&1 | sed 's/^.*staramr //'"), emit: versions_staramr, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gzipped = genome_fasta.getName().endsWith(".gz") ? true : false
    def genome_uncompressed_name = genome_fasta.getName().replace(".gz", "")
    """
    if [ "$is_gzipped" = "true" ]; then
        gzip -c -d $genome_fasta > $genome_uncompressed_name
    fi

    staramr \\
        search \\
        $args \\
        --nprocs $task.cpus \\
        -o ${prefix}_results \\
        $genome_uncompressed_name
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_results
    touch ${prefix}_results/results.xlsx
    echo "" | gzip > ${prefix}_results/{summary,detailed_summary,resfinder,pointfinder,plasmidfinder,mlst}.tsv.gz
    echo "" | gzip > ${prefix}_results/settings.txt.gz
    """
}
