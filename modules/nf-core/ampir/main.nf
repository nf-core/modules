process AMPIR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ampir:1.1.0':
        'biocontainers/r-ampir:1.1.0' }"

    input:
    tuple val(meta), path(faa)
    val model
    val min_length
    val min_probability

    output:
    tuple val(meta), path("*.faa"), emit: amps_faa
    tuple val(meta), path("*.tsv"), emit: amps_tsv
    tuple val("${task.process}"), val('ampir'), eval('Rscript -e "library(ampir); cat(paste((unlist(packageVersion(\'ampir\'))), collapse = \'.\'))"'), emit: versions_ampir, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$faa" == "${prefix}.faa") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    Rscript - <<'EOF'
    library(ampir)

    input_seqs <- read_faa('${faa}')
    prediction <- predict_amps(input_seqs,${min_length},model = '${model}')
    prediction <- prediction[which(prediction\$prob_AMP >= as.numeric(${min_probability})), ]
    output_seqs <- input_seqs[row.names(prediction), ]
    write.table(prediction, file = "${prefix}.tsv", row.names = FALSE, sep = "\t", quote = FALSE, dec = '.')
    df_to_faa(output_seqs, "${prefix}.faa")
    EOF
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("${faa}" == "${prefix}.faa") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    Rscript - <<'EOF'
    library(ampir)

    t <- file("${prefix}.tsv", "w")
    close(t)

    a <- file("${prefix}.faa", "w")
    close(a)
    EOF
    """
}
