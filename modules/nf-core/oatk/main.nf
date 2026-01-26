process OATK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/oatk:1.0':
        'biocontainers/oatk:1.0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(mito_hmm_files)
    tuple val(meta3), path(pltd_hmm_files)

    output:
    tuple val(meta), path("*mito.ctg.fasta"), emit: mito_fasta, optional: true
    tuple val(meta), path("*pltd.ctg.fasta"), emit: pltd_fasta, optional: true
    tuple val(meta), path("*mito.ctg.bed")  , emit: mito_bed, optional: true
    tuple val(meta), path("*pltd.ctg.bed")  , emit: pltd_bed, optional: true
    tuple val(meta), path("*mito.gfa")      , emit: mito_gfa, optional: true
    tuple val(meta), path("*pltd.gfa")      , emit: pltd_gfa, optional: true
    tuple val(meta), path("*annot_mito.txt"), emit: annot_mito_txt, optional: true
    tuple val(meta), path("*annot_pltd.txt"), emit: annot_pltd_txt, optional: true
    tuple val(meta), path("*utg.final.gfa") , emit: final_gfa, optional: true
    tuple val(meta), path("*utg.gfa")       , emit: initial_gfa, optional: true
    tuple val(meta), path("*.log")          , emit: log
    tuple val("${task.process}"), val('oatk'), eval("oatk --version 2>&1"), topic: versions, emit: versions_oatk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mito_hmm_arg = mito_hmm_files ? '-m ' + mito_hmm_files.find { hmm -> hmm.getExtension() =~ /fam|hmm/ } : ""
    def pltd_hmm_arg = pltd_hmm_files ? '-p ' + pltd_hmm_files.find { hmm -> hmm.getExtension() =~ /fam|hmm/ } : ""
    """
    oatk \\
        $args \\
        $mito_hmm_arg \\
        $pltd_hmm_arg \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        ${reads} \\
        2>| >(tee ${prefix}.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annot_mito.txt
    touch ${prefix}.annot_pltd.txt

    touch ${prefix}.mito.gfa
    touch ${prefix}.mito.bed
    touch ${prefix}.mito.ctg.bed
    touch ${prefix}.mito.ctg.fasta

    touch ${prefix}.pltd.gfa
    touch ${prefix}.pltd.bed
    touch ${prefix}.pltd.ctg.bed
    touch ${prefix}.pltd.ctg.fasta

    touch ${prefix}.utg.gfa
    touch ${prefix}.utg.final.gfa

    touch ${prefix}.log
    """
}
