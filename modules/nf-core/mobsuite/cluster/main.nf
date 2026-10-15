process MOBSUITE_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta), path(mobtyper_report), path(taxonomy)
    tuple val(meta2), path(ref_cluster), path(ref_fasta)

    output:
    tuple val(meta), path("${prefix}/clusters.txt"), emit: clusters
    tuple val(meta), path("${prefix}/references_updated.fasta"), emit: fasta
    tuple val(meta), path("${prefix}/references_updated.fasta.msh"), emit: mash
    tuple val(meta), path("${prefix}/references_updated.fasta.n*"), emit: blastdb
    tuple val("${task.process}"), val('mobsuite'), eval("mob_cluster --version | sed 's/mob_cluster //'"), topic: versions, emit: versions_mobsuite

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mode = ref_cluster ? 'update' : 'build'
    def ref_cluster_cmd = ref_cluster ? "--ref_cluster_file ${ref_cluster}" : ''
    def ref_fasta_cmd = ref_fasta ? "--ref_fasta_file ${ref_fasta}" : ''
    """
    mob_cluster \\
        --mode ${mode} \\
        --infile ${fasta} \\
        --mob_typer_file ${mobtyper_report} \\
        --taxonomy ${taxonomy} \\
        ${ref_cluster_cmd} \\
        ${ref_fasta_cmd} \\
        --num_threads ${task.cpus} \\
        --outdir ${prefix} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}
    touch ${prefix}/clusters.txt
    touch ${prefix}/references_updated.fasta
    touch ${prefix}/references_updated.fasta.msh
    touch ${prefix}/references_updated.fasta.ndb
    touch ${prefix}/references_updated.fasta.nhr
    touch ${prefix}/references_updated.fasta.nin
    touch ${prefix}/references_updated.fasta.njs
    touch ${prefix}/references_updated.fasta.not
    touch ${prefix}/references_updated.fasta.nsq
    touch ${prefix}/references_updated.fasta.ntf
    touch ${prefix}/references_updated.fasta.nto
    """
}
