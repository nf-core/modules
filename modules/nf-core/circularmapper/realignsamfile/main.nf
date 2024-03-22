process CIRCULARMAPPER_REALIGNSAMFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circularmapper:1.93.5--h2a3209d_3':
        'biocontainers/circularmapper:1.93.5--h2a3209d_3' }"

    input:
    tuple val(meta), path(reference)
    val(elong)

    output:
    tuple val(meta), path("*_${elong}.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
# default, no filtering. Input SAM does not have to be sorted or anything else, output will be BAM already, saving you to convert ;-)
java -jar RealignSAMFile.jar -e 500 -i output_circmapper.sam -r reference.fasta
# with filtering, same command
java -jar RealignSAMFile.jar -e 500 -i output_circmapper.sam -r reference.fasta -f true -x false

#This will create a file called "output_circmapper_realigned.bam" in the same folder than your input automatically.
samtools sort -@ 4 output_circmapper_realigned.BAM -o output_circmapper_realigned.sorted.bam
samtools index output_circmapper_realigned.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circulargenerator: \$(circulargenerator -h | grep 'usage' | sed 's/usage: CircularGenerator//')
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${elong}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circulargenerator: \$(circulargenerator -h | grep 'usage' | sed 's/usage: CircularGenerator//')
    END_VERSIONS
    """
}
