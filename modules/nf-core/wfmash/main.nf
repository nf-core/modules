process WFMASH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wfmash:0.13.0--h11f254b_0':
        'quay.io/biocontainers/wfmash:0.13.0--h11f254b_0' }"

    input:
    tuple val(meta), path(fasta_gz), path(paf), path(gzi), path(fai)
    val(query_self)
    path(fasta_query_list)

    output:
    tuple val(meta), path("*.paf"), emit: paf
    tuple val("${task.process}"), val('wfmash'), eval('wfmash --version 2>&1 | cut -f 1 -d "-" | cut -f 2 -d "v"'), emit: versions_wfmash, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? task.ext.prefix : paf ? "${meta.id}" + "." + paf.baseName.split("\\.")[-1] : "${meta.id}"
    def query_list = fasta_query_list ? "--query-file-list ${fasta_query_list}" : ""
    def query = query_self ? "${fasta_gz}" : ""
    def paf_mappings = paf ? "--input-paf ${paf}" : ""
    """
    wfmash \\
        ${fasta_gz} \\
        $query \\
        $query_list \\
        --threads $task.cpus \\
        $paf_mappings \\
        $args > ${prefix}.paf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paf
    """
}
