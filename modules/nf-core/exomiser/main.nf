process EXOMISER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::exomiser-rest-prioritiser=13.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/exomiser-rest-prioritiser:13.2.0--hdfd78af_0':
        'biocontainers/exomiser-rest-prioritiser:13.2.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(config)

    # TODO @abhayr20: add additional inputs according to the pattern above
    # See:
    # https://nf-co.re/docs/contributing/modules#new-module-guidelines-and-pr-review-checklist and
    # https://www.youtube.com/watch?v=84XtbqRkKSk&list=PL3xpfTVZLcNikun1FrSvtXW8ic32TciTJ

    output:
    tuple val(meta), path("*.vcf") , optional:true, emit: vcf
    tuple val(meta), path("*.html"), optional:true, emit: html
    tuple val(meta), path("*.json"), optional:true, emit: json
    # TODO @abhayr20: add additional outputs
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #TODO @abhayr20: Add command used to annotate the vcf file
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
