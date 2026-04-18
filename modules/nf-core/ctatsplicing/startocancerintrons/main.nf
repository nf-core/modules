process CTATSPLICING_STARTOCANCERINTRONS {
    tag "$meta.id"
    label 'process_single'

    container "nf-core/ctatsplicing:0.0.3"

    input:
    tuple val(meta), path(split_junction), path(junction), path(bam), path(bai)
    tuple val(meta2), path(genome_lib)

    output:
    tuple val(meta), path("*.cancer_intron_reads.sorted.bam")    , emit: cancer_introns_sorted_bam, optional: true
    tuple val(meta), path("*.cancer_intron_reads.sorted.bam.bai"), emit: cancer_introns_sorted_bai, optional: true
    tuple val(meta), path("*.gene_reads.sorted.sifted.bam")      , emit: gene_reads_sorted_bam    , optional: true
    tuple val(meta), path("*.gene_reads.sorted.sifted.bam.bai")  , emit: gene_reads_sorted_bai    , optional: true
    tuple val(meta), path("*.cancer.introns")                    , emit: cancer_introns
    tuple val(meta), path("*.cancer.introns.prelim")             , emit: cancer_introns_prelim
    tuple val(meta), path("*${prefix}.introns")                  , emit: introns
    tuple val(meta), path("*.introns.for_IGV.bed")               , emit: introns_igv_bed          , optional: true
    tuple val(meta), path("*.ctat-splicing.igv.html")            , emit: igv_html                 , optional: true
    tuple val(meta), path("*.igv.tracks")                        , emit: igv_tracks               , optional: true
    tuple val(meta), path("*.chckpts")                           , emit: chckpts                  , optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('ctatsplicing'), val("0.0.3"), emit: versions_ctatsplicing, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "CTATSPLICING_STARTOCANCERINTRONS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bam_arg = bam ? "--bam_file ${bam}" : ""
    def create_index = bam && !bai ? "samtools index ${bam}" : ""
    """
    ${create_index}

    /usr/local/src/CTAT-SPLICING/STAR_to_cancer_introns.py \\
        --SJ_tab_file ${split_junction} \\
        --chimJ_file ${junction} \\
        ${bam_arg} \\
        --output_prefix ${prefix} \\
        --ctat_genome_lib ${genome_lib} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def create_igv_files = args.contains("--vis") ? "touch ${prefix}.introns.for_IGV.bed && touch ${prefix}.ctat-splicing.igv.html && touch ${prefix}.igv.tracks" : ""
    """
    ${create_igv_files}
    touch ${prefix}.cancer_intron_reads.sorted.bam
    touch ${prefix}.cancer_intron_reads.sorted.bam.bai
    touch ${prefix}.gene_reads.sorted.sifted.bam
    touch ${prefix}.gene_reads.sorted.sifted.bam.bai
    touch ${prefix}.cancer.introns
    touch ${prefix}.cancer.introns.prelim
    touch ${prefix}.introns
    touch ${prefix}.chckpts
    """
}
