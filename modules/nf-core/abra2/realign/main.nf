
process ABRA2_REALIGN {
    tag "$meta.id" + (binding.variables["meta2"] ? ",$meta2.id" : ",none") + (binding.variables["meta3.id"] ? ",$meta3.id" : ",none")
    label 'process_high'

    conda (params.enable_conda ? "bioconda::abra2=2.24" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quay.io/biocontainers/abra2:2.24--h9f5acd7_1':
        'quay.io/biocontainers/abra2:2.24--h9f5acd7_1' }"

    input:
    tuple val(meta), path(bam),path(bam_bai)
    tuple val(meta2), path(bam2),path(bam_bai2)
    tuple val(meta3), path(bam3),path(bam_bai3)
    tuple path(fasta), path(fasta_fai)
    path(targets)

    output:
    tuple val(meta), path("${meta.id}.abra.bam"), emit: bam
    tuple val(meta2), path("${meta2.id}.abra.bam"), optional : true, emit: bam2
    tuple val(meta3), path("${meta3.id}.abra.bam"), optional : true, emit: bam3
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '2.24' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def targets_opt = targets ? "--targets ${targets[0]}" : ""

    def bam_in = bam.getName()
    def bam2_in = meta2.id ? ','+bam2.getName() : ''
    def bam3_in = meta3.id ? ','+bam3.getName() : ''


    def bam_out = meta.id + ".abra.bam"
    def bam2_out = meta2 ? ','+meta2.id+".abra.bam" : ''
    def bam3_out = meta3 ? ','+meta3.id+".abra.bam" : ''

    def single_end = meta.single_end ? '--single' : ''

    """
    abra2 \\
        --threads ${task.cpus} \\
        --tmpdir . \\
        --in ${bam_in}${bam2_in}${bam3_in} \\
        --out ${bam_out}${bam2_out}${bam3_out} \\
        --ref ${fasta} \\
        ${single_end} \\
        ${targets_opt} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abra2: $VERSION
    END_VERSIONS
    """
}
