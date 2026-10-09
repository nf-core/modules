// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLASMIDID {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/plasmidid:1.6.5--hdfd78af_0'
        : 'quay.io/biocontainers/plasmidid:1.6.5--hdfd78af_0'}"

    input:
    tuple val(meta), path(scaffold)
    path fasta

    output:
    tuple val(meta), path("${prefix}/*final_results.html"), emit: html
    tuple val(meta), path("${prefix}/*final_results.tab"), emit: tab
    tuple val(meta), path("${prefix}/images/"), emit: images
    tuple val(meta), path("${prefix}/logs/"), emit: logs
    tuple val(meta), path("${prefix}/data/"), emit: data
    tuple val(meta), path("${prefix}/database/"), emit: database
    tuple val(meta), path("${prefix}/fasta_files/"), emit: fasta_files
    tuple val(meta), path("${prefix}/kmer/"), emit: kmer
    tuple val("${task.process}"), val('plasmidid'), eval("plasmidID --version 2>&1 | sed 's/^plasmidID //'"), topic: versions, emit: versions_plasmidid
    tuple val("${task.process}"), val('parallel'), eval("parallel --version | sed '1!d;s/GNU parallel //'"), topic: versions, emit: versions_parallel
    tuple val("${task.process}"), val('blast'), eval("blastn -version | sed '2!d;s/.*blast //;s/,.*//'"), topic: versions, emit: versions_blast
    tuple val("${task.process}"), val('prokka'), eval("prokka --version 2>&1 | sed 's/prokka //'"), topic: versions, emit: versions_prokka
    tuple val("${task.process}"), val('circos'), eval("circos --version | sed 's/.*v //;s/ .*//'"), topic: versions, emit: versions_circos
    tuple val("${task.process}"), val('mash'), eval("mash --version"), topic: versions, emit: versions_mash

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plasmidID \\
        -d ${fasta} \\
        -s ${prefix} \\
        -c ${scaffold} \\
        ${args} \\
        -o .

    mv NO_GROUP/${prefix} ./${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/images ${prefix}/logs ${prefix}/data ${prefix}/database ${prefix}/fasta_files ${prefix}/kmer
    touch ${prefix}/${prefix}_final_results.html
    touch ${prefix}/${prefix}_final_results.tab
    touch ${prefix}/images/${prefix}_convergence.png
    touch ${prefix}/logs/plasmidID.log
    touch ${prefix}/data/${prefix}.bed
    touch ${prefix}/database/${prefix}_database.fasta
    touch ${prefix}/fasta_files/${prefix}_convergence.fasta
    touch ${prefix}/kmer/${prefix}_kmer.txt
    """
}
