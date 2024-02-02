process ANGSD_GL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/angsd:0.940--hce60e53_2':
        'biocontainers/angsd:0.940--hce60e53_2' }"

    input:
    tuple val(meta), path(bam)
    path(fasta) //Optionally

    output:
    tuple val(meta), path("${prefix}"), emit: genotype_likelihood
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def method = args.contains("-GL 1") ?: 1 : args.contains("-GL 2") : 2 : args.contains("-GL 3") : 3 : args.contains("-GL 4") : 4
    def minq = (method = 3) ?: 0 : 13 //Required for having minq set on 0 for SOAPsnp

    """
    ls -1 *.bam > bamlist.txt

    //GL 1 + 2 work differently than 3 and 4
    if [[${GL} != 3 || ${GL} != 4]]; then
        angsd \\
            -nThreads ${task.cpus} \\
            $args \\
            -GL ${method} \\
            -bam bamlist.txt \\
            -out ${prefix} \\
            -ref ${fasta} \\
            -minq $minq
    fi

    //SOAPsnp needs to run twice (above, then this part again), so checking if GL = 3 requesting SOAP SNP model
    if [${GL} == 3]; then
        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            -GL ${method} \\
            -doGlf 1 \\
            -out ${prefix}
    fi

    //SYK model
        if [${GL} == 4]; then
        angsd \\
            -nThreads ${task.cpus} \\
            -bam bamlist.txt \\
            -GL ${method} \\
            -doGlf 1 \\
            -doCounts 1 \\
            -errors ${prefix}.error \\
            -out ${prefix} \\
            $minq
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ls -1 *.bam > bamlist.txt
    touch ${prefix}
    touch ${prefix}.error

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        angsd: \$(echo \$(angsd 2>&1) | grep version | head -n 1 | sed 's/.*version: //g;s/ .*//g')
    END_VERSIONS
    """
}
