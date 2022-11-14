
process ONCOCNV {
    tag "$tumor_dataset_id"
    label 'process_medium'

    // conda package not yet available
    conda (params.enable_conda ? "bioconda::oncocnv=7.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/oncocnv/v7.0_cv1/oncocnv_v7.0_cv1.sif':
        'registry.hub.docker.com/biocontainers/oncocnv:v7.0_cv1' }"

    input:
    tuple normal_dataset_id, file(normal_bams), file(normal_bais)
    tuple tumor_dataset_id, file(tumor_bams), file(tumor_bais)
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
    def prefix = params.enable_conda ? '' : '/usr/local/bin'
    def cghseg = task.ext.cghseg ? 'cghseg' : ''
    """
    perl ${prefix}/ONCOCNV_getCounts.pl \\
        getControlStats \\
        -m ${task.mode} \\
        -b ${bed} \\
        -c ${normal_bams.join(",")} \\
        -o ControlStats.txt

    perl ${prefix}/ONCOCNV_getCounts.pl \\
        getSampleStats \\
        -m ${task.mode} \\
        -c ControlStats.txt \\
        -s ${tumor_bams.join(",")} \\
        -o SampleStats.txt
    
    cat ControlStats.txt \\
        | grep -v start \\
        | awk '{print $1,$2,$3}' \\
        | sed "s/ /\t/g" > target.bed

    perl ${prefix}/createTargetGC.pl \\
        -bed target.bed \\
        -fi ${fasta} \\
        -od . \\
        -of TargetGC.txt

    cat ${prefix}/processControl.R \\
        | R \\
        --slave \\
        --args ControlStats.txt ControlStatsProcessed.txt target.GC.txt

    cat ${prefix}/processSamples.R \\
        | R \\
        --slave \\
        --args TestStats.txt ControlStatsProcessed.txt Output.log ${cghseg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncocnv: 7.0
        perl: \$(perl --version | grep 'This is perl' | sed 's/.*(//g' | sed 's/)//g')
        r: \$(R --version | grep "R version")
    END_VERSIONS
    """
}
