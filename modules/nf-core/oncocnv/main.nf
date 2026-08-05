process ONCOCNV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1f/1f04dcb708961e66aade70ebd0f9faa1a2980e2eff66bdbd8c947575c8ba4f31/data':
        'community.wave.seqera.io/library/oncocnv_r-pscbs:ffe800b4a6974d56' }"

    input:
    tuple val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index)
    path bed
    path fasta

    output:
    tuple val(meta), path("*.profile.png"), emit: png
    tuple val(meta), path("*.profile.txt"), emit: profile
    tuple val(meta), path("*.summary.txt"), emit: summary
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('oncocnv'), val("7.0"), topic: versions, emit: versions_oncocnv
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed '2!d;s/.*v//;s/).*//'"), topic: versions, emit: versions_perl
    tuple val("${task.process}"), val('r-base'), eval("R --version | sed '1!d; s/.*version //; s/ .*//'"), topic: versions, emit: versions_r_base

    when:
    task.ext.when == null || task.ext.when

    script:
    def cghseg = task.ext.args2 ?: 'cghseg'
    def mode = task.ext.args ?: '-m Ampli'
    def normal_id = normal.join(',')
    def tumor_id = tumor.join(',')

    """
    export R_USER_CACHE_DIR=\$(pwd)/R_cache

    perl \$(which ONCOCNV_getCounts.pl) \\
        getControlStats \\
        ${mode} \\
        -b ${bed} \\
        -c ${normal_id} \\
        -o ControlStats.txt

    perl \$(which ONCOCNV_getCounts.pl) \\
        getSampleStats \\
        ${mode} \\
        -c ControlStats.txt \\
        -s ${tumor_id} \\
        -o SampleStats.txt

    cat ControlStats.txt \\
        | grep -v start \\
        | awk '{print \$1,\$2,\$3}' \\
        | sed "s/ /\t/g" > target.bed

    perl \$(which createTargetGC.pl) \\
        -bed target.bed \\
        -fi ${fasta} \\
        -od . \\
        -of TargetGC.txt

    cat \$(which processControl.R) \\
        | R \\
        --slave \\
        --args ControlStats.txt ControlStatsProcessed.txt TargetGC.txt

    cat \$(which processSamples.R) \\
        | R \\
        --slave \\
        --args SampleStats.txt ControlStatsProcessed.txt Output.log ${cghseg}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile.png
    touch ${prefix}.profile.txt
    touch ${prefix}.summary.txt
    """
}
