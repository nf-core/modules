// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process TRYCYCLER_SUBSAMPLE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trycycler:0.5.3--pyhdfd78af_0':
        'biocontainers/trycycler:0.5.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*/*.fastq.gz") , emit: subreads
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    trycycler \\
        subsample \\
        $args \\
        --reads ${reads} \\
        --threads $task.cpus \\
        --out_dir ${prefix}

    gzip $args2 ${prefix}/*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trycycler: \$(echo \$(trycycler --version 2>&1) | sed 's/Trycycler v\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    count=\$(echo $args | sed -n 's/.*--count\\s*\\([0-9]*\\).*/\\1/p')

    if [[ -z "${count}" ]]; then count=12; fi

    for n in \$(seq -f "%02g" 1 $count); do
        mkdir -p ${prefix}
        touch ${prefix}/sample_${n}.fastq.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trycycler: \$(echo \$(trycycler --version 2>&1) | sed 's/Trycycler v\$//' ))
    END_VERSIONS
    """
}
