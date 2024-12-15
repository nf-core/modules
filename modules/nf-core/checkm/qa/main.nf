process CHECKM_QA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm-genome:1.2.3--pyhdfd78af_1' :
        'biocontainers/checkm-genome:1.2.3--pyhdfd78af_1' }"

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
    def coverage = coverage_file && coverage_file.isFile()             ? "--coverage_file ${coverage_file}"  : ""
    def exclude  = exclude_marker_file && exclude_marker_file.isFile() ? "--exclude_markers ${exclude_marker_file}" : ""
    """
    checkm \\
        qa \\
        --threads ${task.cpus} \\
        --file ${prefix}.${suffix} \\
        ${marker_file} \\
        ${analysis_dir} \\
        ${coverage} \\
        ${exclude} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
