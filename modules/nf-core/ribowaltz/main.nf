// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

process RIBOWALTZ {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribowaltz:2.0--r43hdfd78af_0':
        'biocontainers/ribowaltz:2.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf), path(fasta)

    output:
    tuple val(meta), path("*.psite_offset.tsv.gz"), emit: offset, optional: true
    tuple val(meta), path("offset_plot/*"), emit: offset_plot, optional: true
    tuple val(meta), path("*.psite.tsv.gz"), emit: psites, optional: true
    tuple val(meta), path("*.codon_coverage_rpf.tsv.gz"), emit: codon_coverage_rpf, optional: true
    tuple val(meta), path("*.codon_coverage_psite.tsv.gz"), emit: codon_coverage_psite, optional: true
    tuple val(meta), path("*.cds_coverage_psite.tsv.gz"), emit: cds_coverage, optional: true
    tuple val(meta), path("*nt_coverage_psite.tsv.gz"), emit: cds_window_coverage, optional: true
    tuple val(meta), path("ribowaltz_qc/*.pdf"), emit: ribowaltz_qc, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    template 'ribowaltz.r'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extensions = ['length_bins_for_psite.pdf', 'ends_heatmap.pdf', 'frames.pdf', 'metaprofile_psite.pdf', 'psite_region.pdf'] // Subset of key plots extensions
    def VERSION = '2.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.psite_offset.tsv.gz ${prefix}.psite.tsv.gz ${prefix}.codon_coverage_rpf.tsv.gz ${prefix}.codon_coverage_psite.tsv.gz ${prefix}.cds_coverage_psite.tsv.gz
    mkdir -p offset_plot/${prefix} && touch offset_plot/${prefix}/28.pdf
    mkdir -p ribowaltz_qc 
    for ext in "\${extensions[@]}"; do touch ribowaltz_qc/${prefix}.\${ext}; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        riboWaltz: $VERSION
    END_VERSIONS
    """
}
