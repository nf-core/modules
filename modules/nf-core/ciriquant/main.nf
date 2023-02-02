process CIRIQUANT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::ciriquant=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ciriquant:1.1.2--pyhdfd78af_2' :
        'quay.io/biocontainers/ciriquant:1.1.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(reads)
    path gtf
    path fasta
    path bwa
    path hisat2

    output:
    tuple val(meta), path("${prefix}/${prefix}.gtf"), emit: gtf
    path("${prefix}")                               , optional: true, emit: intermediates_directory
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    bwa_prefix = fasta.toString() == 'genome.fa' ? fasta.toString() : fasta.toString() - ~/.(fa|fasta)$/
    hisat2_prefix = fasta.toString() - ~/.(fa|fasta)$/
    fasta_path = fasta.toRealPath()
    gtf_path = gtf.toRealPath()
    bwa_path = bwa.toRealPath()
    hisat2_path = hisat2.toRealPath()
    """
    BWA=`which bwa`
    HISAT2=`which hisat2`
    STRINGTIE=`which stringtie`
    SAMTOOLS=`which samtools`

    touch travis.yml
    printf "name: ciriquant\ntools:\n  bwa: \$BWA\n  hisat2: \$HISAT2\n  stringtie: \$STRINGTIE\n  samtools: \$SAMTOOLS\n\nreference:\n  fasta: ${fasta_path}\n  gtf: ${gtf_path}\n  bwa_index: ${bwa_path}/${bwa_prefix}\n  hisat_index: ${hisat2_path}/${hisat2_prefix}" >> travis.yml

    CIRIquant \\
        -t ${task.cpus} \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --config travis.yml \\
        -o ${prefix} \\
        -p ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ciriquant: \$(echo \$(CIRIquant --version 2>&1) | sed 's/CIRIquant //g' )
    END_VERSIONS
    """
}
