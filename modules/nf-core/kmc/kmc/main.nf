process KMC_KMC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kmc=3.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmc:3.2.1--hf1761c0_2':
        'biocontainers/kmc:3.2.1--hf1761c0_2' }"

    input:
    tuple val(meta) , path(reads)

    output:
    path "versions.yml"           , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory.toGiga() ?: '2'
    def kmer_size = task.ext.kmer_size
    def input_formats = 0
    def format_str =""
    
    // def extensions = reads.collect { filename -> filename.split(".")[0] }
    // println extensions

    def input_str = reads.collect{filename -> filename.toString()}.join('\n')
    
    """
    echo $input_str > input.txt

    kmc \\
        -sf$task.cpus \\
        -sp$task.cpus \\
        -sr$task.cpus \\
        -t$task.cpus \\
        -m$memory \\
        -k${kmer_size}  \\
        $format_str \\
        $args  \\
        @input.txt \\
        ${kmer_size}.kmc\\
        out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc : \$(echo \$(kmc --version 2>&1 | grep "K-Mer " | awk '{print \$5}' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc : \$(echo \$(kmc --version  | grep "K-Mer " ))
    END_VERSIONS
    """
}
