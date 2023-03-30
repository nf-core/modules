process WFMASH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::wfmash=0.10.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wfmash:0.10.2--hfdddef0_0':
        'quay.io/biocontainers/wfmash:0.10.2--hfdddef0_0' }"

    input:
    tuple val(meta), path(fasta_gz)
    val(query_self)
    path(gzi)
    path(fai)
    path(fasta_query_list)
    path(paf)

    output:
    tuple val(meta), path("*.paf"), emit: paf
    path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when

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


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wfmash: \$(echo \$(wfmash --version 2>&1) | cut -f 1 -d '-' | cut -f 2 -d 'v')
    END_VERSIONS
    """
}
