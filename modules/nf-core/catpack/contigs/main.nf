process CATPACK_CONTIGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_1'
        : 'biocontainers/cat:6.0.1--hdfd78af_1'}"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta2), path(database)
    tuple val(meta3), path(taxonomy)
    tuple val(meta4), path(proteins)
    tuple val(meta5), path(diamond_table)

    output:
    tuple val(meta), path("*.ORF2LCA.txt"), emit: orf2lca
    tuple val(meta), path("*.contig2classification.txt"), emit: contig2classification
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.diamond"), optional: true, emit: diamond
    tuple val(meta), path("*.predicted_proteins.faa"), optional: true, emit: faa
    tuple val(meta), path("*.gff"), optional: true, emit: gff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def premade_proteins = proteins ? "--proteins_fasta ${proteins}" : ''
    def premade_table = diamond_table ? "--diamond_alignment ${diamond_table}" : ''
    """
    CAT_pack contigs \\
        --nproc ${task.cpus} \\
        --contigs_fasta ${contigs} \\
        --database_folder ${database} \\
        --taxonomy_folder ${taxonomy} \\
        --out_prefix ${prefix} \\
        ${premade_proteins} \\
        ${premade_table} \\
        ${args}

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
