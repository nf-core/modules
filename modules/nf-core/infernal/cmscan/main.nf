process INFERNAL_CMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56cfd2ffd3ea93fb3611175e7383ca253108ad32e868bc08067bf9197da408af/data' :
        'community.wave.seqera.io/library/infernal:1.1.5--8877beae3740ff72' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(cmdb)
    path clanin

    output:
    tuple val(meta), path("*.txt.gz"), emit: txt
    tuple val(meta), path("*.tbl.gz"), emit: tbl
    tuple val("${task.process}"), val('infernal'), eval("cmsearch -h | sed '2!d;s/.*INFERNAL //;s/ .*//'"), emit: versions_infernal, topic: versions
    tuple val("${task.process}"), val('gzip'), eval("gzip --version |& sed '1!d;s/gzip //'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"

    // Infernal requirement: --clanin must accompany --fmt 2 and --tblout
    def clan_opts   = clanin ? "--clanin ${clanin} --fmt 2 --oclan" : ""
    def tblout_arg  = "--tblout ${prefix}.tbl"

    """
    cmscan \\
        ${args} \\
        --cpu ${task.cpus} \\
        -o ${prefix}.txt \\
        ${tblout_arg} \\
        ${clan_opts} \\
        ${cmdb}/*.cm \\
        ${fasta}

    gzip --no-name ${prefix}.txt ${prefix}.tbl
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt ${prefix}.tbl
    gzip --no-name ${prefix}.txt ${prefix}.tbl
    """
}
