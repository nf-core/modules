process TXIMETA_TXIMPORT {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta%3A1.20.1--r43hdfd78af_0' :
        'biocontainers/bioconductor-tximeta:1.20.1--r43hdfd78af_0' }"

    input:
    tuple val(meta), path("quants/*")
    tuple val(meta2), path(tx2gene)
    val quant_type

    output:
    tuple val(meta), path("*gene*.tsv")             , emit: gene
    tuple val(meta), path("*transcript*.tsv")       , emit: transcript
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tximport.r'

    stub:
    """
    touch ${meta.id}.gene_tpm.tsv
    touch ${meta.id}.gene_counts.tsv
    touch ${meta.id}.gene_counts_length_scaled.tsv
    touch ${meta.id}.gene_counts_scaled.tsv
    touch ${meta.id}.gene_lengths.tsv
    touch ${meta.id}.transcript_tpm.tsv
    touch ${meta.id}.transcript_counts.tsv
    touch ${meta.id}.transcript_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
