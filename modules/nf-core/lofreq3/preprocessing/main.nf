
process LOFREQ3_PREPROCESSING {
    tag "$meta.id"
    label 'process_low'

    //conda "lofreq:3.0"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/lofreq:3.0':
    //    'quay.io/biocontainers/lofreq:3.0' }"
    container 'quay.io/vojalu/lofreq:3.0'


    input:
    tuple val(meta),
          path(reffa),
          path(fq1),
          path(fq2)

    output:
    tuple val(meta), 
          path("*.bam"), 
          path(".bai"), 
          path(reffa),            , emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    fname = fname = fq1.name.replace('_1','').replace('_R1','').replace('_r1','').replace('.fq','').replace('.fastq','').replace('.gz','')
    obam = fname + ".bam"


    """
    bwa index $reffa

    bwa mem $reffa $fq1 $fq2 | \
        samtools fixmate - - | \
        lofreq viterbi -f $reffa -b - | \
        samtools sort - | \
        lofreq indelqual -f $reffa -b - | \
        lofreq alnqual -f $reffa  -b - | \
        samtools view -b - -o $obam;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq3: \$(echo "3.0 - reimplementation of 2 in Nim")
    END_VERSIONS
    """
}
