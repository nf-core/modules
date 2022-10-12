process MELT {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::melt=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/melt:1.0.3':
        'quay.io/biocontainers/melt:1.0.3' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(transposon_file)
    path(genes_file)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    java -Xmx8G -jar MELT.jar Single \\
        -b hs37d5/NC_007605 \\
        -t ${transposon_file}  \\
        -h ${fasta} \\
        -bamfile $prefix \\
        -w $prefix \\
        -n ${genes_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        melt: \$(echo \$(melt --version 2>&1) | sed 's/^.*melt //; s/ .*\$//' ))
    END_VERSIONS
    """
}
