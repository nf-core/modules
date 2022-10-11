process SMOOVE {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::smoove=0.2.8" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/smoove:0.2.8--h9ee0642_0':
        'quay.io/biocontainers/smoove:0.2.8--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    smoove call \\
	--outdir . \\
	--name ${prefix} \\
	--fasta ${fasta} \\
	-p $task.cpus \\
	${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smoove: \$(echo \$(smoove --version 2>&1) | sed 's/^.*smoove //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
