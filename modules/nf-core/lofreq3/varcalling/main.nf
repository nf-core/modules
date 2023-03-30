
process LOFREQ3_VARCALLING {
    tag "$meta.id"
    label 'process_low'

    //conda "lofreq:3.0"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/lofreq:3.0':
    //    'quay.io/biocontainers/lofreq:3.0' }"
    container 'quay.io/vojalu/lofreq:3.0'
"

    input:
    tuple val(meta),
          path(reffa),
          path(bam),
          path(bai),
    val pileup

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (pileup) pileup_command = "-p"
    else pileup_command = ""

    """
    lofreq call $pileup_command -b $bam -f $reffa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         lofreq3: \$(echo "3.0 - reimplementation of 2 in Nim")
    END_VERSIONS
    """
}
