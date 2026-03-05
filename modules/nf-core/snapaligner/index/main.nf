process SNAPALIGNER_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap-aligner:2.0.5--h077b44d_2':
        'biocontainers/snap-aligner:2.0.5--h077b44d_2' }"

    input:
    tuple val(meta), path(fasta), path(altcontigfile), path(nonaltcontigfile), path(altliftoverfile)

    output:
    tuple val(meta), path("snap/*") ,emit: index
    tuple val("${task.process}"), val('snap-aligner'), eval("snap-aligner 2>&1 | sed 's/^.*version //;s/.\$//;q'"), topic: versions, emit: versions_snapaligner

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def altcontigfile_arg = altcontigfile ? '-altContigFile ' + altcontigfile : ''
    def nonaltcontigfile_arg = nonaltcontigfile ? '-nonAltContigFile ' + nonaltcontigfile : ''
    def altliftoverfile_arg = altliftoverfile ? '-altLiftoverFile ' + altliftoverfile : ''
    """
    mkdir snap

    snap-aligner \\
        index \\
        $fasta \\
        snap \\
        -t${task.cpus} \\
        $altcontigfile_arg \\
        $nonaltcontigfile_arg \\
        $altliftoverfile_arg \\
        $args
    """
    stub:
    """
    mkdir snap
    echo "Genome" > snap/Genome
    echo "GenomeIndex" > snap/GenomeIndex
    echo "GenomeIndexHash" > snap/GenomeIndexHash
    echo "OverflowTable" > snap/OverflowTable
    """
}
