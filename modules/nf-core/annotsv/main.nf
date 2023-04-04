process ANNOTSV {
    tag "$meta.id"
    label 'process_low'

    // Creation of the bioconda recipe is underway, but awaits some vital changes to the source code (https://github.com/lgmgeo/AnnotSV/issues/166)
    container 'quay.io/cmgg/annotsv:3.3.2'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "ANNOTSV module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple val(meta), path(variants), path(variants_index)
    tuple val(meta2), path(annotations)
    tuple val(meta3), path(candidate_genes)
    tuple val(meta4), path(candidate_small_variants)
    tuple val(meta5), path(false_positive_snv)
    tuple val(meta6), path(gene_transcripts)

    output:
    tuple val(meta), path("*.tsv")              , emit: tsv
    tuple val(meta), path("*.unannotated.tsv")  , emit: unannotated_tsv, optional:true
    // VCF output is a bit flaky and is prone to failures (https://github.com/lgmgeo/AnnotSV/issues/168)
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
        -outputFile ${prefix}.tsv \\
        -SVinputFile ${variants} \\
        ${args}

    mv *_AnnotSV/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotsv: \$(echo \$(AnnotSV -help) 2>&1 | head -n1 | sed 's/^AnnotSV //')
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
        annotsv: \$(echo \$(AnnotSV -help) 2>&1 | head -n1 | sed 's/^AnnotSV //')
    END_VERSIONS
    """
}
