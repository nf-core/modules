// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process MOBSUITE_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7ca9972cb410cbd9e028437e0dd286bb76986565610432ecbb8e38ac63744bc4/data'
        : 'community.wave.seqera.io/library/mob_suite:3.1.9--eaef3d58d6b4b3f2'}"

    input:
    tuple val(meta), path(fasta), path(mobtyper_report), path(taxonomy)
    tuple val(meta2), path(ref_cluster, stageAs: 'ref/*'), path(ref_fasta, stageAs: 'ref/*')
    path(taxa_sqlite, stageAs: 'taxa/*')

    output:
    tuple val(meta), path("${prefix}/clusters.txt"), emit: clusters
    tuple val(meta), path("${prefix}/references_updated.fasta"), emit: fasta
    tuple val(meta), path("${prefix}/references_updated.fasta.msh"), emit: mash
    tuple val(meta), path("${prefix}/references_updated.fasta.n*"), emit: blastdb
    tuple val("${task.process}"), val('mobsuite'), eval("mob_cluster --version | sed 's/mob_cluster //'"), topic: versions, emit: versions_mobsuite

    when:
    task.ext.when == null || task.ext.when

    script:
    if (!ref_cluster != !ref_fasta) {
        error("ref_cluster and ref_fasta must be provided together (update mode) or both omitted (build mode)")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def mode = ref_cluster ? 'update' : 'build'
    def ref_cluster_cmd = ref_cluster ? "--ref_cluster_file ${ref_cluster}" : ''
    def ref_fasta_cmd = ref_fasta ? "--ref_fasta_file ${ref_fasta}" : ''
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    mkdir -p db
    ln -s ../${taxa_sqlite} db/taxa.sqlite

    mob_cluster \\
        --mode ${mode} \\
        --infile ${fasta_name} \\
        --mob_typer_file ${mobtyper_report} \\
        --taxonomy ${taxonomy} \\
        ${ref_cluster_cmd} \\
        ${ref_fasta_cmd} \\
        --database_directory db \\
        --num_threads ${task.cpus} \\
        --outdir ${prefix} \\
        ${args}

    [ -s ${prefix}/clusters.txt ] || { echo "mob_cluster failed, see log above" >&2; exit 1; }
    """

    stub:
    if (!ref_cluster != !ref_fasta) {
        error("ref_cluster and ref_fasta must be provided together (update mode) or both omitted (build mode)")
    }
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
