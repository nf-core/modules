process VG_DECONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::vg=1.38.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.38.0--h9ee0642_0' :
        'quay.io/biocontainers/vg:1.38.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fasta), path(gfa)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg deconstruct \\
        --threads $task.cpus \\
        --path=$fasta \\
        $args \\
        $gfa > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | grep -o 'vg .*' | cut -f3 -d ' ')
    END_VERSIONS
    """
}
