process ANNOTSV_ANNOTSV {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::annotsv=3.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/annotsv:3.3.4--py311hdfd78af_1' :
        'quay.io/biocontainers/annotsv:3.3.4--py311hdfd78af_1' }"

    input:
    tuple val(meta), path(variants), path(variants_index)
    tuple val(meta2), path(annotations)
    tuple val(meta3), path(candidate_genes)
    tuple val(meta4), path(candidate_small_variants)
    tuple val(meta5), path(false_positive_snv)
    tuple val(meta6), path(gene_transcripts)

    output:
    tuple val(meta), path("*.tsv")              , emit: tsv
    tuple val(meta), path("*.unannotated.tsv")  , emit: unannotated_tsv, optional: true
    tuple val(meta), path("*.vcf")              , emit: vcf, optional:true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def cand_genes = candidate_genes ? "-candidateGenesFile ${candidate_genes}" : ""
    def cand_small = candidate_small_variants ? "-candidateSnvIndelFiles ${candidate_small_variants}" : ""
    def fp_snv = false_positive_snv ? "-snvIndelFiles ${false_positive_snv}" : ""
    def transcripts = gene_transcripts ? "-txFile ${gene_transcripts}" : ""

    """
    AnnotSV \\
        -annotationsDir ${annotations} \\
        ${cand_genes} \\
        ${cand_small} \\
        ${fp_snv} \\
        ${transcripts} \\
        -outputFile ${prefix}.tsv \\
        -SVinputFile ${variants} \\
        ${args}

    mv *_AnnotSV/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def create_vcf = args.contains("-vcf 1") ? "touch ${prefix}.vcf" : ""

    """
    touch ${prefix}.tsv
    touch ${prefix}.unannotated.tsv
    ${create_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help 2>&1 | head -n1 | sed 's/^AnnotSV //'))
    END_VERSIONS
    """
}
