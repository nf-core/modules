process ANGSD_SOAPSNPCALIBRATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--h13024bc_4':
        'quay.io/biocontainers/angsd:0.940--h13024bc_4' }"

    input:
    tuple val(meta), path(bams), path(bam_indices)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(meta), path("*_calibration_matrix"), emit: gl_calibration
    tuple val("${task.process}"), val('angsd'), eval("angsd 2>&1 | sed '1!d;s/.*version: //;s/ .*//'"), emit: versions_angsd, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''

    // Touch fai index to ensure it is newer than the fasta (ANGSD requirement)
    def touch_ref = reference_fai  ? "sleep 1 && touch ${reference_fai}"  : ''
    
    """
    ${touch_ref}
    printf '%s\\n' ${bams} > bamlist.txt

    angsd \\
        -nThreads ${task.cpus} \\
        -bam bamlist.txt \\
        -minQ 0 \\
        -GL 3 \\
        -ref ${reference_fasta} \\
        -out ${prefix} \\
        ${args} \\
        -tmpdir ${prefix}_calibration_matrix
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_calibration_matrix
    touch ${prefix}_calibration_matrix/0.counts0
    touch ${prefix}_calibration_matrix/0.qual0
    touch ${prefix}_calibration_matrix/1.counts1
    touch ${prefix}_calibration_matrix/1.qual1
    touch ${prefix}_calibration_matrix/2.counts2
    touch ${prefix}_calibration_matrix/2.qual2
    touch ${prefix}_calibration_matrix/3.counts3
    touch ${prefix}_calibration_matrix/3.qual3
    """
}
