process YLEAF {
    tag "$meta.id"
    label 'process_medium'

    // TODO AH: Remove this once yleaf 3.3.0 is released to bioconda
    conda "nf-core"

    input:
    tuple val(meta), path(input_file)
    path reference_fasta, stageAs: "reference.fa"
    path y_reference_fasta, stageAs: "y_reference.fa"
    val reference_genome
    val reads_threshold
    val quality_thresh
    val base_majority
    val prediction_quality
    val draw_haplogroups
    val collapsed_draw_mode
    val ancient_dna
    val private_mutations
    val minor_allele_frequency

    output:
    tuple val(meta), path("${meta.id}/hg_prediction.hg"), emit: haplogroup
    tuple val(meta), path("${meta.id}/*.log"), emit: log
    tuple val(meta), path("${meta.id}/hg_tree_image.pdf"), emit: tree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_genome = reference_genome ?: "hg38"
    // These defaults correspond to the defaults used in the Yleaf package
    def reads_thresh = reads_threshold ?: 10
    def qual_thresh = quality_thresh ?: 20
    def base_maj = base_majority ?: 90
    def pred_qual = prediction_quality ?: 0.95
    def draw_hg = draw_haplogroups ?: false
    def collapsed_mode = collapsed_draw_mode ?: false
    def ancient = ancient_dna ?: false
    def private_mut = private_mutations ?: false
    def maf = minor_allele_frequency ?: 0.01

    """
    mkdir -p ${prefix}

    echo "Input file: ${input_file}"
    echo "Reference genome: ${ref_genome}"

    PYTHONPATH=\$PWD:/home/a/Yleaf Yleaf \\
        --vcffile ${input_file} \\
        --output ${prefix} \\
        --reference_genome ${ref_genome} \\
        --full_genome_reference ${reference_fasta} \\
        --y_chromosome_reference ${y_reference_fasta} \\
        --threads ${task.cpus} \\
        --force \\
        --reanalyze \\
        --reads_treshold ${reads_thresh} \\
        --quality_thresh ${qual_thresh} \\
        --base_majority ${base_maj} \\
        --prediction_quality ${pred_qual} \\
        ${draw_hg ? '--draw_haplogroups' : ''} \\
        ${collapsed_mode ? '--collapsed_draw_mode' : ''} \\
        ${ancient ? '--ancient_DNA' : ''} \\
        ${private_mut ? '--private_mutations' : ''} \\
        --minor_allele_frequency ${maf} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: \$(Yleaf --help | grep "version" | sed 's/.*version //g' | sed 's/).*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    touch ${prefix}/run.log
    touch ${prefix}/hg_prediction.hg
    touch ${prefix}/hg_tree_image.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: "3.3.0"
    END_VERSIONS
    """
}
