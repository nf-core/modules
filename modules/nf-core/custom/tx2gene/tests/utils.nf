process REARRANGE {

    input:
    tuple val(meta), path ("quants")

    output:
    tuple val(meta), path("quants/*", includeInputs: true), emit: subdirs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    """
}
