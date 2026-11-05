// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TCOFFEE_ALNCOMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e85c972da17c3e7d5e8b283cf61b106111180e3fc02d4b04f0933106b4ca50b3/data'
        : 'community.wave.seqera.io/library/t-coffee_pigz:aa6b4f1ff4adb194'}"

    input:
    tuple val(meta), path(msa), path(ref_msa)

    output:
    tuple val(meta), path("*.scores"), emit: scores
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ? task.ext.prefix : "${msa.baseName}"
    def args = task.ext.args ? task.ext.args.contains('compare_mode') ? task.ext.args  : (task.ext.args +'-compare_mode tc ' ) : '-compare_mode tc '
    def metric_name = args.split('compare_mode ')[1].split(' ')[0]
    def header = meta.keySet().join(",")
    def values = meta.values().join(",")

    """
    # check whether it is compressed
    if [[ "${msa}" == *.gz ]]; then
        unpigz -c ${msa} > uncompressed_msa.fa
    else
        mv ${msa} uncompressed_msa.fa
    fi

    export TEMP='./'
    t_coffee -other_pg aln_compare \
        -al1 ${ref_msa} \
        -al2 uncompressed_msa.fa \
        ${args} \
        | grep -v "seq1" | grep -v '*' | \
        awk '{ print \$4}' ORS="\t" \
        >> "scores.txt"

    # Add metadata info to output file
    echo "${header},${metric_name}" > "${prefix}.scores"

    # Add values
    scores=\$(awk '{sub(/[[:space:]]+\$/, "")} 1' scores.txt | tr -s '[:blank:]' ',')
    echo "${values},\$scores" >> "${prefix}.scores"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch "${prefix}.scores"
    """
}
