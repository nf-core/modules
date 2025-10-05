process DIANN_PRELIMINARYANALYSIS {
    tag "$ms_file.baseName"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/diann/v1.8.1_cv1/diann_v1.8.1_cv1.img' :
        'docker.io/biocontainers/diann:v1.8.1_cv1' }"

    input:
    tuple val(meta), path(ms_file), path(predict_library)

    output:
    tuple val(meta), path("*.quant"), emit: diann_quant
    tuple val(meta), path("*_diann.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN_PRELIMINARYANALYSIS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''

    """
    diann   --lib ${predict_library} \\
            --f ${ms_file} \\
            --threads ${task.cpus} \\
            --temp ./ \\
            ${args}

    cp report.log.txt ${ms_file.baseName}_diann.log

    # Check for mzML indexing errors first (specific format issue)
    if grep -q -E "(No index list offset found|Failure reading input file)" ${ms_file.baseName}_diann.log; then
        echo "ERROR: Invalid mzML format - missing index offset or corrupted file"
        echo "DIA-NN requires fully indexed mzML files with proper indexListOffset and valid MS2 spectra"
        exit 1
    fi

    # Check for other serious errors
    if grep -q -E "(ERROR.*failed to load.*files:|ERROR.*cannot load the file.*skipping|No MS2 spectra.*aborting)" ${ms_file.baseName}_diann.log; then
        echo "ERROR: DIA-NN PRELIMINARY ANALYSIS failed - check log for file loading errors"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DIA-NN: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "DIANN_PRELIMINARYANALYSIS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''

    """
    touch ${ms_file.baseName}.quant
    touch ${ms_file.baseName}_diann.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DIA-NN: \$(diann 2>&1 | grep "DIA-NN" | grep -oP "\\d+\\.\\d+(\\.\\w+)*(\\.[\\d]+)?")
    END_VERSIONS
    """
}
