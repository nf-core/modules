process ALLELECOUNTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cancerit-allelecount=4.3.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cancerit-allelecount:4.3.0--h41abebc_0' :
        'quay.io/biocontainers/cancerit-allelecount:4.3.0--h41abebc_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    path loci
    path fasta

    output:
    tuple val(meta), path("*.alleleCount"), emit: allelecount
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_options = fasta ? "-r $fasta": ""

    """
    alleleCounter \\
        $args \\
        -l $loci \\
        -b $input \\
        $reference_options \\
        -o ${prefix}.alleleCount

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        allelecounter: \$(alleleCounter --version)
    END_VERSIONS
    """
}
