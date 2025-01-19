process CATPACK_CONTIGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_0'
        : 'biocontainers/cat:6.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta2), path(database)
    tuple val(meta3), path(taxonomy)

    output:
    tuple val(meta), path("*.ORF2LCA.txt"), emit: orf2lca
    tuple val(meta), path("*.contig2classification.txt"), emit: contig2classification
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.diamond"), emit: diamond
    tuple val(meta), path("*.predicted_proteins.faa"), emit: faa
    tuple val(meta), path("*.gff"), emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    CAT_pack contigs \\
        -n ${task.cpus} \\
        -c ${contigs} \\
        -d ${database} \\
        -t ${taxonomy} \\
        -o ${prefix}
        ${args} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ORF2LCA.txt
    touch ${prefix}.contig2classification.txt
    touch ${prefix}.log
    touch ${prefix}.diamond
    touch ${prefix}.predicted_proteins.faa
    touch ${prefix}.predicted_proteins.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """
}
