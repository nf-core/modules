process MULTIQC_SAV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d745674ad0e37f8eccbf52c0e2c1f0e9aaebc9dd4b548bca803db08e800b7d3/data' :
        'community.wave.seqera.io/library/multiqc_multiqc_sav:b7474a8ef88d6535' }"

    input:
    tuple val(meta), path(runinfo_xml), path(interop_bin, stageAs: "InterOp/*")
    path(extra_multiqc_files, stageAs: "?/*")
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    path "*.html"      , emit: report
    path "*_data"      , emit: data
    path "*_plots"     , optional:true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions_multiqc
    tuple val("${task.process}"), val('multiqc_sav'), eval('python -c "import multiqc_sav; print(multiqc_sav.__version__)"'), emit: versions_multiqc_sav
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    def extra_config = extra_multiqc_config ? "--config ${extra_multiqc_config}" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${extra_config} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}
