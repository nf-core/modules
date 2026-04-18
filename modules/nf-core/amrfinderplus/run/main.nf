process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.7--hf69ffd2_0':
        'biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    env 'VER'                                       , emit: tool_version
    env 'DBVER'                                     , emit: db_version
    tuple val("${task.process}"), val("amrfinderplus"), eval("amrfinder --version"), emit: versions_amrfinderplus, topic: versions
    tuple val("${task.process}"), val("amrfinderplus_database"), eval("amrfinder --database amrfinderdb --database_version 2>&1 | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}\\.[0-9]+' | tail -1"), emit: versions_amrfinderplus_database, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def is_compressed_fasta = fasta.getName().endsWith(".gz") ? true : false
    def is_compressed_db = db.getName().endsWith(".gz") ? true : false

    organism_param = meta.containsKey("organism") ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fasta_name = fasta.getName().replace(".gz", "")
    fasta_param = "-n"
    if (meta.containsKey("is_proteins")) {
        if (meta.is_proteins) {
            fasta_param = "-p"
        }
    }
    """
    if [ "${is_compressed_fasta}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if [ "${is_compressed_db}" == "true" ]; then
        mkdir amrfinderdb
        tar xzvf ${db} -C amrfinderdb
    else
        mv ${db} amrfinderdb
    fi

    amrfinder \\
        ${fasta_param} ${fasta_name} \\
        ${organism_param} \\
        ${args} \\
        --database amrfinderdb \\
        --threads ${task.cpus} > ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(amrfinder --database amrfinderdb --database_version 2>&1 | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}\\.[0-9]+' | tail -1)
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: DB version information will not be provided for the database in the stub run, since it is not actually downloaded through the update module
    """
    echo ${args}

    touch ${prefix}.tsv
    ${meta.containsKey("organism") ? "touch ${prefix}-mutations.tsv" : ""}

    VER=\$(amrfinder --version)
    DBVER="2026-01-01.1"
    """
}
