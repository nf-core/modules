process RIBOWALTZ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:2.0--r43hdfd78af_0':
        'biocontainers/ribowaltz:2.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.best_offset.txt")                , emit: best_offset          , optional: true
    tuple val(meta), path("*.psite_offset.tsv{,.gz}")         , emit: offset               , optional: true
    tuple val(meta), path("offset_plot/*")                    , emit: offset_plot          , optional: true
    tuple val(meta), path("*.psite.tsv{,.gz}")                , emit: psites               , optional: true
    tuple val(meta), path("*.codon_coverage_rpf.tsv{,.gz}")   , emit: codon_coverage_rpf   , optional: true
    tuple val(meta), path("*.codon_coverage_psite.tsv{,.gz}") , emit: codon_coverage_psite , optional: true
    tuple val(meta), path("*.cds_coverage_psite.tsv{,.gz}")   , emit: cds_coverage         , optional: true
    tuple val(meta), path("*nt_coverage_psite.tsv{,.gz}")     , emit: cds_window_coverage  , optional: true
    tuple val(meta), path("ribowaltz_qc/*.pdf")               , emit: ribowaltz_qc         , optional: true
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'ribowaltz.r'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.best_offset.txt
    touch ${prefix}.psite_offset.tsv
    touch ${prefix}.psite.tsv
    touch ${prefix}.codon_coverage_rpf.tsv
    touch ${prefix}.codon_coverage_psite.tsv
    touch ${prefix}.cds_coverage_psite.tsv
    mkdir -p offset_plot/${prefix} && touch offset_plot/${prefix}/29.pdf
    mkdir -p ribowaltz_qc && touch ribowaltz_qc/${prefix}.metaprofile_psite.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-ribowaltz: \$(Rscript -e "library(riboWaltz); cat(as.character(packageVersion('riboWaltz')))")
    END_VERSIONS
    """
}
