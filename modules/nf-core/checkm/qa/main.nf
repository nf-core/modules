process CHECKM_QA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::checkm-genome=1.2.1" : null)
    def container_image = "checkm-genome:1.2.1--pyhdfd78af_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(analysis_dir), path(marker_file), path(coverage_file)
    path exclude_marker_file

    output:
    tuple val(meta), path("${prefix}.txt")  , optional: true, emit: output
    tuple val(meta), path("${prefix}.fasta"), optional: true, emit: fasta
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    suffix       = task.ext.args?.matches(".*-o 9.*|.*--out_file 9.*") ? "fasta" : "txt"
    def coverage = coverage_file                                   ? "--coverage_file ${coverage_file}"  : ""
    def exclude  = exclude_marker_file                             ? "--exclude_markers ${marker_filer}" : ""
    """
    checkm \\
        qa \\
        --threads ${task.cpus} \\
        --file ${prefix}.${suffix} \\
        $marker_file \\
        $analysis_dir \\
        $coverage \\
        $exclude \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
