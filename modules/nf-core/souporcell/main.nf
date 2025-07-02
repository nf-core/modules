process SOUPORCELL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker.io/nicotru/souporcell:0.0.4"

    input:
    tuple val(meta), path(bam), path(barcodes)
    tuple val(meta2), path(fasta)
    val(clusters)

    output:
    tuple val(meta), path("test/*.vcf"), emit: vcf
    tuple val(meta), path("test/*.tsv"), emit: tsv
    path "versions.yml", emit: versions

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
        souporcell: \$(souporcell_pipeline.py --version 2>&1 || echo "N/A")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    touch ${prefix}/fake.vcf
    touch ${prefix}/fake.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell: \$(souporcell --version)
    END_VERSIONS
    """
}
