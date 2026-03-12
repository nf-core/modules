process METABULI_BUILD {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.1.1--pl5321h0bb26bb_0':
        'biocontainers/metabuli:1.1.1--pl5321h0bb26bb_0' }"

    input:
    tuple val(meta), path(fasta)
    path taxonomy_names, stageAs: 'taxonomy/names.dmp'
    path taxonomy_nodes, stageAs: 'taxonomy/nodes.dmp'
    path taxonomy_merged, stageAs: 'taxonomy/merged.dmp'
    path accession2taxid, stageAs: 'taxonomy/*'
    path cds_info

    output:
    tuple val(meta), path("$prefix"), emit: db
    tuple val("${task.process}"), val('metabuli'), eval('metabuli 2>&1 | awk \'/metabuli Version:/ {print $3}\''), emit: versions_metabuli, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    make_merged = taxonomy_merged ? "" : "touch taxonomy/merged.dmp"
    cds_info_arg = cds_info ? "--cds-info cds_info.txt" : ""
    """
    $make_merged
    echo $fasta | tr ' ' '\\n' > fasta.txt
    echo $cds_info | tr ' ' '\\n' > cds_info.txt

    metabuli build \\
        "${prefix}" \\
        fasta.txt \\
        $accession2taxid \\
        --taxonomy-path taxonomy \\
        --max-ram ${task.memory.toGiga()} \\
        --threads ${task.cpus} \\
        ${cds_info_arg} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "$prefix"

    touch "$prefix/acc2taxid.map"
    touch "$prefix/diffIdx"
    touch "$prefix/info"
    touch "$prefix/split"
    touch "$prefix/taxID_list"
    touch "$prefix/db.parameters"
    """
}
