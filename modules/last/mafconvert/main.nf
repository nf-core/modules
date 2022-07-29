process LAST_MAFCONVERT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::last=1250' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/last:1250--h2e03b76_0' :
        'quay.io/biocontainers/last:1250--h2e03b76_0' }"

    input:
    tuple val(meta), path(maf)
    val(format)

    output:
    tuple val(meta), path("*.axt.gz"),      optional:true, emit: axt_gz
    tuple val(meta), path("*.blast.gz"),    optional:true, emit: blast_gz
    tuple val(meta), path("*.blasttab.gz"), optional:true, emit: blasttab_gz
    tuple val(meta), path("*.chain.gz"),    optional:true, emit: chain_gz
    tuple val(meta), path("*.gff.gz"),      optional:true, emit: gff_gz
    tuple val(meta), path("*.html.gz"),     optional:true, emit: html_gz
    tuple val(meta), path("*.psl.gz"),      optional:true, emit: psl_gz
    tuple val(meta), path("*.sam.gz"),      optional:true, emit: sam_gz
    tuple val(meta), path("*.tab.gz"),      optional:true, emit: tab_gz
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    maf-convert $args $format $maf | gzip --no-name \\
        > ${prefix}.${format}.gz

    # maf-convert has no --version option but lastdb (part of the same package) has.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        last: \$(lastdb --version 2>&1 | sed 's/lastdb //')
    END_VERSIONS
    """
}
