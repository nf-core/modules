process SAM2LCA_UPDATEDB {
    tag "${acc2tax_name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.4--pyhdfd78af_0':
        'biocontainers/sam2lca:1.1.4--pyhdfd78af_0' }"

    input:
    val(acc2tax_name)
    val(taxo_db_name)
    path(taxo_nodes)   // nodes.dmp
    path(taxo_names)   // names.dmp
    path(taxo_merged)  // merged.dmp
    path(acc2tax_json) // optional
    path(acc2tax)      // acc2tax.gz
    path(acc2tax_md5)  // acc2tax.gz.md5

    output:
    path (prefix)       , emit: sam2lca_db
    path "versions.yml" , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    prefix      = task.ext.prefix ?: 'sam2lca_db'
    def args    = task.ext.args ?: ''
    def names   = taxo_names ? "--taxo_names ${taxo_names}" : ''
    def nodes   = taxo_nodes ? "--taxo_nodes ${taxo_nodes}" : ''
    def merged  = taxo_merged ? "--taxo_merged ${taxo_merged}" : ''
    def json    = acc2tax_json ? "--acc2tax_json ${acc2tax_json}" : ''
    """
    mkdir -p $prefix

    sam2lca -d $prefix \\
        update-db \\
        -t $taxo_db_name \\
        $names \\
        $nodes \\
        $merged \\
        -a $acc2tax_name \\
        $json \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/^sam2lca, version //' )
    END_VERSIONS
    """

    stub:
    prefix = 'sam2lca_db'
    """
    mkdir -p $prefix
    touch ${prefix}/test.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/^sam2lca, version //' )
    END_VERSIONS
    """
}
