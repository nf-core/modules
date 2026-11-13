process EUKCC_SINGLE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/40/4007a76cef091bbb1d959ad608d94fb5bac8fa7aea96bca9df83cb72445555de/data':
        'quay.io/biocontainers/eukcc:2.3.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}/eukcc.tsv")                , emit: tsv
    tuple val(meta), path("${prefix}/eukcc.log")                , emit: log
    tuple val(meta), path("${prefix}/scmg_marker_table.csv.gz") , emit: marker_table, optional: true
    tuple val("${task.process}"), val('eukcc'), eval("eukcc --version | sed 's/^EukCC version //'"), topic: versions, emit: versions_eukcc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    eukcc \\
        single \\
        ${args} \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --out ${prefix} \\
        ${fasta}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    printf "fasta\\tcompleteness\\tcontamination\\tncbi_lng\\tncbi_lng_names\\n" > ${prefix}/eukcc.tsv
    touch ${prefix}/eukcc.log
    """
}
