process VARNET_FILTER {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varnet:1.5.4--pyhdfd78af_0':
        'quay.io/biocontainers/varnet:1.5.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_tumor), path(index_tumor), path(input_normal), path(index_normal)
    tuple val(meta2), path(intervals)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)
    tuple val(meta5), path(resource_dir)

    output:
    tuple val(meta), path("${prefix}/candidates"), emit: candidates
    tuple val("${task.process}"), val("varnet"), eval("varnet-filter --version | sed 's/VarNet //'"), emit: versions_varnet, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def regions = intervals ? "--region_bed ${intervals}" : ""
    def normal = input_normal ? "--normal_bam ${input_normal}" : ""
    def resources = resource_dir ? "--resource_dir ${resource_dir}" : ""
    if (!input_normal && !resource_dir) {
        error "VARNET_FILTER tumor-only mode needs the germline filtering resources folder (download it with varnet-download-resources, see https://doi.org/10.6084/m9.figshare.34420128)."
    }
    """
    export TF_CPP_MIN_LOG_LEVEL=3

    varnet-filter \\
        --sample_name ${prefix} \\
        ${normal} \\
        --tumor_bam ${input_tumor} \\
        --reference ${fasta} \\
        --output_dir . \\
        --processes ${task.cpus} \\
        ${regions} \\
        ${resources} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/candidates/snvs ${prefix}/candidates/indels
    touch ${prefix}/candidates/snvs/Positions.csv
    touch ${prefix}/candidates/indels/Positions.csv
    """
}
