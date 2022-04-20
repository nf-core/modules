process BOWTIE2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4 bioconda::samtools=1.15.1 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0' :
        'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:1744f68fe955578c63054b55309e05b41c37a80d-0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    val   save_unaligned

    output:
    tuple val(meta), path('*.bam')    , emit: bam
    tuple val(meta), path('*.log')    , emit: log
    tuple val(meta), path('*fastq.gz'), emit: fastq, optional:true
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        def unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
        bowtie2 \\
            -x \$INDEX \\
            -U $reads \\
            --threads $task.cpus \\
            $unaligned \\
            $args \\
            2> ${prefix}.bowtie2.log \\
            | samtools view -@ $task.cpus $args2 -bhS -o ${prefix}.bam -

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
        END_VERSIONS
        """
    } else {
        def unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
        bowtie2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --threads $task.cpus \\
            $unaligned \\
            $args \\
            2> ${prefix}.bowtie2.log \\
            | samtools view -@ $task.cpus $args2 -bhS -o ${prefix}.bam -

        if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
            mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
        fi
        if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
            mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
        END_VERSIONS
        """
    }
}
