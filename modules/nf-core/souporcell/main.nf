process SOUPORCELL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92c054bdfc9170bd58c09de480160923d35dd67008650733f3d03588520082b1/data' :
        'community.wave.seqera.io/library/souporcell:2.5--2b23aea4d0753391' }"

    input:
    tuple val(meta), path(bam), path(barcodes), val(clusters)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*/clusters.tsv")         , emit: clusters
    tuple val(meta), path("*/cluster_genotypes.vcf"), emit: vcf
    tuple val(meta), path("*/ambient_rna.txt")      , emit: ambient_rna
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    """
    mkdir -p temp
    export TMPDIR=./temp
    souporcell_pipeline.py \\
        -i $bam \\
        -b $barcodes \\
        -f $fasta \\
        -t $task.cpus \\
        -o $prefix \\
        -k $clusters \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell: 2.5
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    touch ${prefix}/clusters.tsv
    touch ${prefix}/cluster_genotypes.vcf
    touch ${prefix}/ambient_rna.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell: 2.5
    END_VERSIONS
    """
}
