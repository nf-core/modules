process STRINGTIE_MERGE {
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda     (params.enable_conda ? "bioconda::stringtie=2.2.1" : null)
    def container_image = "stringtie:2.2.1--hecb563c_2"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    path stringtie_gtf
    path annotation_gtf

    output:
    path "stringtie.merged.gtf", emit: gtf
    path  "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    """
    stringtie \\
        --merge $stringtie_gtf \\
        $reference \\
        -o stringtie.merged.gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch stringtie.merged.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
