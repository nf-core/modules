process ONCOCNV {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::oncocnv=7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/oncocnv/v7.0_cv2/oncocnv_v7.0_cv2.sif':
        'registry.hub.docker.com/biocontainers/oncocnv:v7.0_cv2' }"

    input:
    tuple val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index)
    path bed
    path fasta

    output:
    tuple val(meta), path("*.profile.png")  ,emit: png
    tuple val(meta), path("*.profile.txt")  ,emit: profile
    tuple val(meta), path("*.summary.txt")  ,emit: summary
    path "versions.yml"                     ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cghseg = task.ext.args2 ?: 'cghseg'
    def mode = task.ext.args ?: '-m Ampli'
    def normal = normal.join(',')
    def tumor = tumor.join(',')
    def VERSION = '7.0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    perl \$(which ONCOCNV_getCounts.pl) \\
        getControlStats \\
        $mode \\
        -b ${bed} \\
        -c $normal \\
        -o ControlStats.txt

    perl \$(which ONCOCNV_getCounts.pl) \\
        getSampleStats \\
        $mode \\
        -c ControlStats.txt \\
        -s $tumor \\
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
        --args SampleStats.txt ControlStatsProcessed.txt Output.log $cghseg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncocnv: $VERSION
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(v//g' | sed 's/)//g')
        r: \$(R --version | grep "R version" | sed 's/R version //g')
    END_VERSIONS
    """
}
