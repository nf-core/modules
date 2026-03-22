process HDL_RG {
    tag "hdl_rg_${meta.id}_${meta2.id}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-base_r-data.table_r-hdl:a2553504418194cc' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-hdl:cb9d70356e12d034' }"

    input:
    tuple val(meta), path(sumstats1)
    tuple val(meta2), path(sumstats2)
    tuple val(meta3), path(hdl_ref_panel_dir)

    output:
    tuple val(meta), val(meta2), path("${meta.id}.${meta2.id}.rg.tsv"), emit: correlation_results
    tuple val(meta), val(meta2), path("${meta.id}.${meta2.id}.hdl.log"), emit: rg_log
    tuple val("${task.process}"), val("hdl"), val("1.4.0"), emit: versions_hdl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'hdl_rg.R'

    stub:
    """
    cat <<EOF_STUB > ${meta.id}.${meta2.id}.rg.tsv
trait1\ttrait2\trg\tse\tp\th2_trait1\th2_trait2\tcovariance\teigen_use
${meta.id}\t${meta2.id}\t0.3500\t0.0700\t0.00001\t0.2000\t0.1800\t0.0600\tautomatic
EOF_STUB

    cat <<EOF_STUB > ${meta.id}.${meta2.id}.hdl.log
Stub HDL rg run for ${meta.id} and ${meta2.id}
EOF_STUB
    """
}
