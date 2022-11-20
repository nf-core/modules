
process ONCOCNV {
    tag "$tumor_dataset_id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::oncocnv=7.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/oncocnv/v7.0_cv1/oncocnv_v7.0_cv1.sif':
        'registry.hub.docker.com/biocontainers/oncocnv:v7.0_cv1' }"

    input:
    tuple val(normal_dataset_id), path(normal_bams), path(normal_bais)
    tuple val(tumor_dataset_id), path(tumor_bams), path(tumor_bais)
    path bed
    path fasta

    output:
    tuple path("*.profile.png")  ,emit: png
    tuple path("*.profile.txt")  ,emit: profile
    tuple path("*.summary.txt")  ,emit: summary
    path "versions.yml"          ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cghseg = task.ext.cghseg ? 'cghseg' : ''
    def mode = task.mode ?: 'Ampli'
    def normal_bams_input = normal_bams.join(',')
    def prefix = params.enable_conda ? '/opt/conda/bin' : '/usr/local/bin'
    def tumor_bams_input = tumor_bams.join(',')
    """
    ls /usr/share/miniconda
    ls /usr/share/miniconda/*
    perl ${prefix}/ONCOCNV_getCounts.pl \\
        getControlStats \\
        -m $mode \\
        -b ${bed} \\
        -c $normal_bams_input \\
        -o ControlStats.txt

    perl ${prefix}/ONCOCNV_getCounts.pl \\
        getSampleStats \\
        -m $mode \\
        -c ControlStats.txt \\
        -s $tumor_bams_input \\
        -o SampleStats.txt

    cat ControlStats.txt \\
        | grep -v start \\
        | awk '{print \$1,\$2,\$3}' \\
        | sed "s/ /\t/g" > target.bed

    perl ${prefix}/createTargetGC.pl \\
        -bed target.bed \\
        -fi ${fasta} \\
        -od . \\
        -of TargetGC.txt

    cat ${prefix}/processControl.R \\
        | R \\
        --slave \\
        --args ControlStats.txt ControlStatsProcessed.txt TargetGC.txt

    cat ${prefix}/processSamples.R \\
        | R \\
        --slave \\
        --args SampleStats.txt ControlStatsProcessed.txt Output.log ${cghseg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncocnv: 7.0
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(v//g' | sed 's/)//g')
        r: \$(R --version | grep "R version" | sed 's/R version //g')
    END_VERSIONS
    """
}
