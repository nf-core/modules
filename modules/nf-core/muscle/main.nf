process MUSCLE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/muscle:3.8.1551--h7d875b9_6' :
        'quay.io/biocontainers/muscle:3.8.1551--h7d875b9_6' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.afa") , emit: aligned_fasta, optional: true
    tuple val(meta), path("*.phyi"), emit: phyi         , optional: true
    tuple val(meta), path("*.phys"), emit: phys         , optional: true
    tuple val(meta), path("*.clw") , emit: clustalw     , optional: true
    tuple val(meta), path("*.html"), emit: html         , optional: true
    tuple val(meta), path("*.msf") , emit: msf          , optional: true
    tuple val(meta), path("*.tree"), emit: tree         , optional: true
    tuple val(meta), path("*.log") , emit: log
    tuple val("${task.process}"), val('muscle'), eval("muscle -version | sed 's/MUSCLE v//;s/ by.*//'"), emit: versions_muscle, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_out = args.contains('-fasta')    ? "-fastaout ${prefix}_muscle_msa.afa" : ''
    def clw_out   = args.contains('-clw')      ? "-clwout ${prefix}_muscle_msa.clw"   : ''
    def msf_out   = args.contains('-msf')      ? "-msfout ${prefix}_muscle_msa.msf"   : ''
    def phys_out  = args.contains('-phys')     ? "-physout ${prefix}_muscle_msa.phys" : ''
    def phyi_out  = args.contains('-phyi')     ? "-phyiout ${prefix}_muscle_msa.phyi" : ''
    def html_out  = args.contains('-html')     ? "-htmlout ${prefix}_muscle_msa.html" : ''
    def tree_out  = args.contains('-maketree') ? "-out ${prefix}_muscle_msa.tree"     : ''
    """
    muscle \\
        ${args} \\
        -in ${fasta} \\
        ${fasta_out} \\
        ${clw_out} \\
        ${msf_out} \\
        ${phys_out} \\
        ${phyi_out} \\
        ${html_out} \\
        ${tree_out} \\
        -loga muscle_msa.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_out = args.contains('-fasta')    ? "touch ${prefix}_muscle_msa.afa" : ''
    def clw_out   = args.contains('-clw')      ? "touch ${prefix}_muscle_msa.clw" : ''
    def msf_out   = args.contains('-msf')      ? "touch ${prefix}_muscle_msa.msf" : ''
    def phys_out  = args.contains('-phys')     ? "touch ${prefix}_muscle_msa.phys" : ''
    def phyi_out  = args.contains('-phyi')     ? "touch ${prefix}_muscle_msa.phyi" : ''
    def html_out  = args.contains('-html')     ? "touch ${prefix}_muscle_msa.html" : ''
    def tree_out  = args.contains('-maketree') ? "touch ${prefix}_muscle_msa.tree" : ''
    """
    ${fasta_out}
    ${clw_out}
    ${msf_out}
    ${phys_out}
    ${phyi_out}
    ${html_out}
    ${tree_out}
    touch muscle_msa.log
    """

}
