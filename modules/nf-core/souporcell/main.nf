process SOUPORCELL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/souporcell_gxx:9e4498a90e1d49fc' :
        'community.wave.seqera.io/library/souporcell_gxx:9f443b06d1671e4e' }"

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
