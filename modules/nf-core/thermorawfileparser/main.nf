process THERMORAWFILEPARSER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/thermorawfileparser:1.4.3--ha8f3691_0' :
        'biocontainers/thermorawfileparser:1.4.3--ha8f3691_0' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("*.{mzML,mzML.gz,mgf,mgf.gz,parquet,parquet.gz}"), emit: spectra
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--format 0") || args.contains("-f 0") ? "mgf" :
                args.contains("--format 1") || args.contains("-f 1") ? "mzML" :
                args.contains("--format 2") || args.contains("-f 2") ? "mzML" :
                args.contains("--format 3") || args.contains("-f 3") ? "parquet" :
                "mzML"
    suffix = args.contains("--gzip")? "${suffix}.gz" : "${suffix}"

    """
    ThermoRawFileParser.sh \\
        --input $raw \\
        --output_file ${prefix}.${suffix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        thermorawfileparser: \$(ThermoRawFileParser.sh --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains("--format 0") || args.contains("-f 0") ? "mgf" :
                args.contains("--format 1") || args.contains("-f 1") ? "mzML" :
                args.contains("--format 2") || args.contains("-f 2") ? "mzML" :
                args.contains("--format 3") || args.contains("-f 3") ? "parquet" :
                "mzML"
    suffix = args.contains("--gzip")? "${suffix}.gz" : "${suffix}"

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        thermorawfileparser: \$(ThermoRawFileParser.sh --version)
    END_VERSIONS
    """
}
