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
    def formatSuffixMap = ["0": "mgf", "1": "mzML", "2": "mzML", "3": "parquet"]
    // Please provide the output format via --format, not -f. If --format not specified, the parser defaults to mzML
    suffix = formatSuffixMap.get(task.ext.args?.format, "mzML")

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
    def formatSuffixMap = ["0": "mgf", "1": "mzML", "2": "mzML", "3": "parquet"]
    // Please provide the output format via --format, not -f. If --format not specified, it defaults to mzML
    suffix = formatSuffixMap.get(task.ext.args?.format, "mzML")

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        thermorawfileparser: \$(ThermoRawFileParser.sh --version)
    END_VERSIONS
    """
}
