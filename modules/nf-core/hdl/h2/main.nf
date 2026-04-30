process HDL_H2 {
    tag "hdl_h2_${meta.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-data.table_r-hdl:a2553504418194cc' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-hdl:cb9d70356e12d034' }"

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta2), path(hdl_ref_panel_dir)

    output:
    tuple val(meta), path("${meta.id}.h2.tsv"), emit: heritability_results
    tuple val(meta), path("${meta.id}.hdl.log"), emit: h2_log
    tuple val("${task.process}"), val("hdl"), val("1.4.0"), emit: versions_hdl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'hdl_h2.R'

    stub:
    """
    cat <<EOF_STUB > ${meta.id}.h2.tsv
trait\th2\tse\tp\teigen_use
${meta.id}\t0.2100\t0.0300\t0.00010\tautomatic
EOF_STUB

    cat <<EOF_STUB > ${meta.id}.hdl.log
Stub HDL h2 run for ${meta.id}
EOF_STUB
    """
}
