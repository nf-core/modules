process FASTAVALIDATOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/py_fasta_validator:0.6--py37h595c7a6_0':
        'biocontainers/py_fasta_validator:0.6--py37h595c7a6_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.success.log')  , emit: success_log , optional: true
    tuple val(meta), path('*.error.log')    , emit: error_log   , optional: true
    tuple val("${task.process}"), val('py_fasta_validator'), eval('py_fasta_validator --version | cut -d" " -f3'), emit: versions_py_fasta_validator, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Fastavalidator module does not support Conda. Please use Docker / Singularity instead."
    }
    def deprecation_message = """
    WARNING: This module has been deprecated.

    Reason:
    This module is no longer recommended for use to validate FASTA files as it is not
    maintained by the original developers.
    It is recommended to use fa-lint
    - nf-core/modules/falint

    """
    assert false: deprecation_message

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    py_fasta_validator \\
        -f $fasta \\
        2> "${prefix}.error.log" \\
        || echo "Errors from fasta_validate printed to ${prefix}.error.log"

    if [ \$(cat "${prefix}.error.log" | wc -l) -gt 0 ]; then
        echo "Validation failed..."

        cat \\
            "${prefix}.error.log"
    else
        echo "Validation successful..."

        mv \\
            "${prefix}.error.log" \\
            fasta_validate.stderr

        echo "Validation successful..." \\
            > "${prefix}.success.log"
    fi
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Fastavalidator module does not support Conda. Please use Docker / Singularity instead."
    }
    def deprecation_message = """
    WARNING: This module has been deprecated.

    Reason:
    This module is no longer recommended for use to validate FASTA files as it is not
    maintained by the original developers.
    It is recommended to use fa-lint
    - nf-core/modules/falint

    """
    assert false: deprecation_message

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "Validation successful..." \\
        > "${prefix}.success.log"
    """
}
