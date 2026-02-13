process LSA_COSINE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4d94f159b95315adf8bf54fdc9db88db10a5aef72dca6245dd163b91e9e0437e/data' :
        'community.wave.seqera.io/library/r-base_r-lsa_r-pheatmap_r-optparse_pruned:901156bc11e60b28' }"

    input:
    tuple val(meta), path(expression_matrix)

    output:
    tuple val(meta), path("*_matrix.csv") , emit: matrix
    tuple val(meta), path("*_heatmap.png"), emit: heatmap
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // HEMOS QUITADO EL 'def' DE AQUÍ ABAJO PARA QUE LA PLANTILLA LO VEA
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'cosine.R'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_matrix.csv
    touch ${prefix}_heatmap.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/ (.*//g')
        r-lsa: 0.73.3
        r-pheatmap: 1.0.12
    END_VERSIONS
    """
}
