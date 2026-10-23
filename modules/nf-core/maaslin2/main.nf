process MAASLIN2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'community.wave.seqera.io/library/maaslin2:1.16.0--68568bbb134347d6' :
        'biocontainers/maaslin2:1.16.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(input_data)
    path input_metadata
    val(output_dir)
    val(fixed_effects)
    val(reference)
    val(min_prevalence)
    val(min_abundance)
    val(normalization)

    output:
    tuple val(meta), path("${output_dir}"), emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(Maaslin2)

    # Load input data
    df_input_data <- read.table("${input_data}",
                                header = TRUE,
                                sep = "\t",
                                row.names = 1,
                                stringsAsFactors = FALSE)

    df_input_metadata <- read.table("${input_metadata}",
                                    header = TRUE,
                                    sep = "\t",
                                    row.names = 1,
                                    stringsAsFactors = FALSE)

    # Run Maaslin2 analysis
    fit_data <- Maaslin2(
        input_data = df_input_data,
        input_metadata = df_input_metadata,
        output = "${output_dir}",
        fixed_effects = c(${fixed_effects}),
        reference = c(${reference}),
        min_prevalence = ${min_prevalence},
        min_abundance = ${min_abundance},
        normalization = "${normalization}",
        $args
    )

    # Output versions
    writeLines(c("\\"${task.process}\\":",
                paste0("    maaslin2: ", packageVersion("Maaslin2"))),
                "versions.yml")
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${output_dir}
    touch ${output_dir}/all_results.tsv
    touch ${output_dir}/significant_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        maaslin2: 1.18.0
    END_VERSIONS
    """
}
