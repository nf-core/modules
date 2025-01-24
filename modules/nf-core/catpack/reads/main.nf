process CATPACK_READS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/cat:6.0.1--hdfd78af_0'
        : 'biocontainers/cat:6.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(contigs)
    tuple val(meta3), path(database)
    tuple val(meta4), path(taxonomy)
    tuple val(meta5), path(bins, stageAs: 'bins/')
    val mode
    tuple val(meta6), path(bam_aligned)
    tuple val(meta7), path(bam_unaligned)
    tuple val(meta8), path(contig2classification)
    tuple val(meta9), path(bin2classification)
    tuple val(meta10), path(unclassified2classification)
    tuple val(meta11), path(proteins)
    tuple val(meta12), path(diamond_alignment)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def insert_reads = meta.single_end ? "--read_file1 ${reads}" : "--read_file1 ${reads[0]} --read_file2 ${reads[1]}"
    def insert_bins = bins ? "--bin_fasta bins/" : ''
    def insert_bam_aligned = bam_aligned ? "-bam1 ${bam_aligned}" : ''
    def insert_bam_unaligned = bam_unaligned ? "--alignment_unmapped ${bam_unaligned}" : ''
    def insert_c2c = contig2classification ? "--c2c ${contig2classification}" : ''
    def insert_b2c = bin2classification ? "--b2c ${bin2classification}" : ''
    def insert_u2c = unclassified2classification ? "--u2c ${unclassified2classification}" : ''
    def insert_proteins = proteins ? "--proteins_fasta ${proteins}" : ''
    def insert_diamond_alignment = diamond_alignment ? "--diamond_alignment ${diamond_alignment}" : ''
    """
    CAT_pack reads \\
        ${insert_reads} \\
        ${args} \\
        -d ${database} \\
        -c ${contigs} \\
        -t ${taxonomy} \\
        -m ${mode} \\
        -o ${prefix}.txt \\
        ${insert_bins} \\
        ${insert_bam_aligned} \\
        ${insert_bam_unaligned} \\
        ${insert_c2c} \\
        ${insert_b2c} \\
        ${insert_u2c} \\
        ${insert_proteins} \\
        ${insert_diamond_alignment}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def insert_reads = meta.single_end ? "--read_file1 ${reads}" : "--read_file1 ${reads[0]} --read_file2 ${reads[1]}"
    def insert_bins = bins ? "-b-bin_fasta ${bins}" : ''
    def insert_bam_aligned = bam_aligned ? "-bam1 ${bam_aligned}" : ''
    def insert_bam_unaligned = bam_unaligned ? "--alignment_unmapped ${bam_unaligned}" : ''
    def insert_c2c = contig2classification ? "--c2c ${contig2classification}" : ''
    def insert_b2c = bin2classification ? "--b2c ${bin2classification}" : ''
    def insert_u2c = unclassified2classification ? "--u2c ${unclassified2classification}" : ''
    def insert_proteins = proteins ? "--proteins_fasta ${proteins}" : ''
    def insert_diamond_alignment = diamond_alignment ? "--diamond_alignment ${diamond_alignment}" : ''
    """
    echo "    CAT_pack reads \\
        ${insert_reads} \\
        ${args} \\
        -d ${database} \\
        -c ${contigs} \\
        -t ${taxonomy} \\
        -m ${mode} \\
        -o ${prefix}.txt \\
        ${insert_bins} \\
        ${insert_bam_aligned} \\
        ${insert_bam_unaligned} \\
        ${insert_c2c} \\
        ${insert_b2c} \\
        ${insert_u2c} \\
        ${insert_proteins} \\
        ${insert_diamond_alignment}"

    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        catpack: \$(CAT_pack --version | sed 's/CAT_pack pack v//g;s/ .*//g')
    END_VERSIONS
    """
}
