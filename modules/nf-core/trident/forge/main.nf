process TRIDENT_FORGE {
    tag "${forge_fn ?: ''} ${forge_s ?: ''}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/poseidon-trident:2.2.2.1--hf7d7819_0'
        : 'quay.io/biocontainers/poseidon-trident:2.2.2.1--hf7d7819_0'}"

    input:
    tuple path(archive_dir, stageAs: 'input_archives/*'), val(forge_s), path(forge_fn)

    output:
    tuple val(out_meta), path("output_dir/${prefix}"), emit: forged_package
    tuple val("${task.process}"), val('trident'), eval('trident --version'), emit: versions_trident, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "forged_package"
    def forge_string = forge_s ? "--forgeString ${forge_s}" : ''
    def forge_file = forge_fn ? "--forgeFile ${forge_fn}" : ''
    def argsParts = args.split(/\s+/)
    out_meta = [id: prefix]
    """
    trident forge \
        -d ${archive_dir} \
        ${forge_file} \
        ${forge_string} \
        -o output_dir/${prefix} \
        -n ${prefix} \
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "forged_package"
    def forge_string = forge_s ? "--forgeString ${forge_s}" : ''
    def forge_file = forge_fn ? "--forgeFile ${forge_fn}" : ''
    def argsParts = args.split(/\s+/)
    out_meta = [id: prefix]
    """
    echo trident forge -d ${archive_dir} ${forge_file} ${forge_string} -o output_dir/${prefix} -n ${prefix} ${args}
    echo "meta: ${out_meta}"
    """
}
