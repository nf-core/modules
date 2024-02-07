
process TXIMETA_TXIMPORT {
    label "process_medium"

    conda "bioconda::bioconductor-tximeta=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-tximeta:1.20.1--r41hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-tximeta:1.20.1--r43hdfd78af_0' }"

    input:
    tuple val(meta), path("quants/*")
    tuple val(meta2), path(tx2gene)
    tuple val(meta3), path(coldata)
    val quant_type

    output:
    path "*gene_tpm.tsv"                 , emit: tpm_gene
    path "*gene_counts.tsv"              , emit: counts_gene
    path "*gene_counts_length_scaled.tsv", emit: counts_gene_length_scaled
    path "*gene_counts_scaled.tsv"       , emit: counts_gene_scaled
    path "*gene_lengths.tsv"             , emit: lengths_gene
    path "*transcript_tpm.tsv"           , emit: tpm_transcript
    path "*transcript_counts.tsv"        , emit: counts_transcript
    path "*transcript_lengths.tsv"       , emit: lengths_transcript
    path "versions.yml"                  , emit: versions

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
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
