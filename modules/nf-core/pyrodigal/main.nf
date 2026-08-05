process PYRODIGAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84cfe1298be8ed1a95ce41684798172a5eb47c3aa172b80f778dd39494cc65e9/data':
        'community.wave.seqera.io/library/pyrodigal_pigz:35ca0f5a299e74f1' }"

    input:
    tuple val(meta), path(fasta)
    val(output_format)

    output:
    tuple val(meta), path("*.${output_format}.gz")      , emit: annotations
    tuple val(meta), path("*.fna.gz")                   , emit: fna
    tuple val(meta), path("*.faa.gz")                   , emit: faa
    tuple val(meta), path("*.score.gz")                 , emit: score
    tuple val("${task.process}"), val('pyrodigal'), eval("pyrodigal --version |& sed 's/pyrodigal v//'"), emit: versions_pyrodigal, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pigz -cdf ${fasta} > pigz_fasta.fna

    pyrodigal \\
        -j ${task.cpus} \\
        $args \\
        -i pigz_fasta.fna \\
        -f $output_format \\
        -o "${prefix}.${output_format}" \\
        -d ${prefix}.fna \\
        -a ${prefix}.faa \\
        -s ${prefix}.score

    pigz -nmf ${prefix}*
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.${output_format}.gz
    echo "" | gzip > ${prefix}.fna.gz
    echo "" | gzip > ${prefix}.faa.gz
    echo "" | gzip > ${prefix}.score.gz
    """
}
