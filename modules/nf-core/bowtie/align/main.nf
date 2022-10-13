process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    def image = "mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:676c5bcfe34af6097728fea60fb7ea83f94a4a5f-0"
    conda (params.enable_conda ? 'bioconda::bowtie=1.3.0 bioconda::samtools=1.15.1' : null)
    container { (params.container_registry ?: 'quay.io' ) + image }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.out'), emit: log
    path  "versions.yml"          , emit: versions
    tuple val(meta), path('*fastq.gz'), optional:true, emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = params.save_unaligned ? "--un ${prefix}.unmapped.fastq" : ''
    def endedness = meta.single_end ? "$reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    INDEX=`find -L ./ -name "*.3.ebwt" | sed 's/.3.ebwt//'`
    bowtie \\
        --threads $task.cpus \\
        --sam \\
        -x \$INDEX \\
        -q \\
        $unaligned \\
        $args \\
        $endedness \\
        2> ${prefix}.out \\
        | samtools view $args2 -@ $task.cpus -bS -o ${prefix}.bam -

    if [ -f ${prefix}.unmapped.fastq ]; then
        gzip ${prefix}.unmapped.fastq
    fi
    if [ -f ${prefix}.unmapped_1.fastq ]; then
        gzip ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
