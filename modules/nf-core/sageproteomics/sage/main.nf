process SAGEPROTEOMICS_SAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sage-proteomics:0.14.7--h031d066_0' :
        'quay.io/biocontainers/sage-proteomics:0.14.7--h031d066_0' }"

    input:
    tuple val(meta),  path("*.mzML")
    tuple val(meta2), path(fasta_proteome)
    tuple val(meta3), path(base_config)

    output:
    tuple val(meta), path("results.sage.tsv"),        emit: results_tsv
    tuple val(meta), path("results.json"),            emit: results_json
    tuple val(meta), path("results.sage.pin"),        emit: results_pin
    tuple val("${task.process}"), val('sageproteomics'), eval("sage --version |& sed '1!d ; s/sage //'"), topic: versions, emit: versions_sageproteomics

    //optional outs
    tuple val(meta), path("tmt.tsv"), optional: true, emit: tmt_tsv
    tuple val(meta), path("lfq.tsv"), optional: true, emit: lfq_tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export RAYON_NUM_THREADS=$task.cpus

    sage $base_config \\
        --disable-telemetry-i-dont-want-to-improve-sage \\
        --fasta $fasta_proteome \\
        --write-pin \\
        *.mzML
    """

    stub:
    """
    touch results.json
    touch results.sage.tsv
    touch results.sage.pin
    """
}
