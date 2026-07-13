nextflow.enable.dsl=2

process blasttblastx {
    tag "${query.getBaseName()}"

    input:
    path query
    path db_dir
    
    output:
    path "tblastx_results.txt", emit: result

    script:
    """
    tblastx \\
        -query ${query} \\
        -db blastdb/nt \\
        -evalue 1e-5 \\
        -outfmt 6 \\
        -out tblastx_results.txt
    """
}

workflow {
Channel
    .fromPath(params.query)
    .set { query }

Channel
    .fromPath(params.db)
    .set { db_dir }

    blasttblastx(query, db_dir)
}
