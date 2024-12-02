process PLINK_GENOME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h031d066_5':
        'biocontainers/plink:1.90b6.21--h031d066_5' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.genome"), emit: genome
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plink \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        --genome \\
        $args \\
        --threads $task.cpus \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
    END_VERSIONS
    """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   touch ${prefix}.genome

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       plink: \$(echo \$(plink --version 2>&1) | sed 's/^PLINK v//' | sed 's/..-bit.*//' )
   END_VERSIONS
   """
   
}
