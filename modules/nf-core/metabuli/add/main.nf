process METABULI_ADD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabuli:1.0.5--pl5321h6a68c12_1':
        'biocontainers/metabuli:1.0.5--pl5321h6a68c12_1' }"

    input:
    tuple val(meta), path(fasta)
    path taxonomy_names, stageAs: 'taxonomy/names.dmp'
    path taxonomy_nodes, stageAs: 'taxonomy/nodes.dmp'
    path taxonomy_merged, stageAs: 'taxonomy/merged.dmp'
    path accession2taxid, stageAs: 'taxonomy/*'

    output:
    tuple val(meta), path("$prefix"), emit: db
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    make_merged = taxonomy_merged ? "" : "touch ${prefix}/taxonomy/merged.dmp"

    """
    mkdir -p $prefix
    mv "taxonomy" $prefix
    $make_merged

    realpath $fasta > fasta.txt
    metabuli add-to-library \\
        fasta.txt \\
        $prefix/$accession2taxid \\
        $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli version)
    END_VERSIONS
    """
}
