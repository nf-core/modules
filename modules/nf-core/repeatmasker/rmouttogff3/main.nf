process REPEATMASKER_RMOUTTOGFF3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(out)

    output:
    tuple val(meta), path("*.gff3") , emit: gff3
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    rm_path=\$(dirname \$(realpath \$(which RepeatMasker)))

    PERL5LIB=\$rm_path rmOutToGFF3.pl \\
        $out \\
        > ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmasker: \$(RepeatMasker -v | sed 's/RepeatMasker version //1')
    END_VERSIONS
    """
}
