// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TCOFFEE_TCS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e8/e85c972da17c3e7d5e8b283cf61b106111180e3fc02d4b04f0933106b4ca50b3/data'
        : 'community.wave.seqera.io/library/t-coffee_pigz:aa6b4f1ff4adb194'}"

    input:
    tuple val(meta) , path(msa)
    tuple val(meta2), path(lib)

    output:
    tuple val(meta), path("*.tcs")   , emit: tcs
    tuple val(meta), path("*.scores"), emit: scores
    tuple val("${task.process}"), val('tcoffee'), eval('t_coffee -version | sed -n \'s/.*Version_\\([^ ]*\\).*/\\1/p\''), emit: versions_tcoffee, topic: versions
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/^.*pigz[[:space:]]*//"'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def lib_arg       = lib ? "-lib ${lib}" : ""
    def header        = meta.keySet().join(",")
    def values        = meta.values().join(",")
    def unzipped_name = msa.toString() - '.gz'
    """
    export TEMP='./'
    filename=${msa}
    if [[ \$(basename $msa) == *.gz ]] ; then
        unpigz -f $msa
        filename=${unzipped_name}
    fi

    # Bad hack to circumvent t_coffee bug
    # Issue described already in: https://github.com/cbcrg/tcoffee/issues/3
    # Add an A in front of filename if the file begins with A
    first_letter_filename=\${filename:0:1}
    if [ "\$first_letter_filename" == "A" ]; then input="A"\$filename; else input=\$filename;  fi

    t_coffee -infile \$input \
        -evaluate -output=score_ascii \
        ${lib_arg} \
        -outfile ${prefix}.tcs

    # Add metadata info to output file
    echo "${header},TCS" > "${prefix}.scores"

    # Add values
    scores=\$(grep 'SCORE=' ${prefix}.tcs | cut -d '=' -f 2 )
    echo "${values},\$scores" >> "${prefix}.scores"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch ${prefix}.tcs
    touch ${prefix}.scores
    """
}
