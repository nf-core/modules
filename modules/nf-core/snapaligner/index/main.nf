process SNAPALIGNER_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::snap-aligner=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snap-aligner:2.0.3--hd03093a_0':
        'biocontainers/snap-aligner:2.0.3--hd03093a_0' }"

    input:
    tuple val(meta), path(fasta), path(altcontigfile), path(nonaltcontigfile), path(altliftoverfile)

    output:
    tuple val(meta), path("snap/*") ,emit: index
    path "versions.yml"             ,emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //')
    END_VERSIONS
    """
    stub:
    """
    mkdir snap
    echo "Genome" > snap/Genome
    echo "GenomeIndex" > snap/GenomeIndex
    echo "GenomeIndexHash" > snap/GenomeIndexHash
    echo "OverflowTable" > snap/OverflowTable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapaligner: \$(snap-aligner 2>&1| head -n 1 | sed 's/^.*version //;s/\.\$//')
    END_VERSIONS
    """
}
