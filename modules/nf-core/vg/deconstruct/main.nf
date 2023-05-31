process VG_DECONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vg=1.43.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.43.0--h9ee0642_0' :
        'biocontainers/vg:1.43.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(gfa)
    path(pb)
    path(gbwt)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def snarls = pb ? "--snarls ${pb}" : ""
    def gbwt_arg = gbwt ? "--gbwt ${gwbt}" : ""
    """
    vg deconstruct \\
        --threads $task.cpus \\
        $args \\
        $snarls \\
        $gbwt_arg \\
        $gfa > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(vg version 2>&1 | grep -o 'vg .*' | cut -f3 -d ' ' | cut -f2 -d 'v')
    END_VERSIONS
    """
}
