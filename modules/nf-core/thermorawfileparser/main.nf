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
    tuple val(meta), path("*.{mzML,mfg,parquet}"), emit: spectra
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Please provide the output format via --format, not -f
    if (task.ext.args?.format) {
        suffix = task.ext.args.format switch {
            case "0" -> "mgf"
            case "1" -> "mzML"
            case "2" -> "mzML"
            case "3" -> "parquet"
            default  -> "mzML"
        }
    } else {
        suffix = "mzML"
    }

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
    // Please provide the output format via --format, not -f
    if (task.ext.args?.format) {
        suffix = task.ext.args.format switch {
            case "0" -> "mgf"
            case "1" -> "mzML"
            case "2" -> "mzML"
            case "3" -> "parquet"
            default  -> "mzML"
        }
    } else {
        suffix = "mzML"
    }

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        thermorawfileparser: \$(ThermoRawFileParser.sh --version)
    END_VERSIONS
    """
}
