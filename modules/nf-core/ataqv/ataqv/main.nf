process ATAQV_ATAQV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ataqv:1.3.1--py310ha155cf9_1' :
        'biocontainers/ataqv:1.3.1--py310ha155cf9_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peak_file)
    val organism
    val mito_name
    path tss_file
    path excl_regs_file
    path autosom_ref_file

    output:
    tuple val(meta), path("*.ataqv.json"), emit: json
    tuple val(meta), path("*.problems")  , emit: problems, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mito = mito_name ? "--mitochondrial-reference-name ${mito_name}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def peak        = peak_file        ? "--peak-file $peak_file"                       : ''
    def tss         = tss_file         ? "--tss-file $tss_file"                         : ''
    def excl_regs   = excl_regs_file   ? "--excluded-region-file $excl_regs_file"       : ''
    def autosom_ref = autosom_ref_file ? "--autosomal-reference-file $autosom_ref_file" : ''
    """
    ataqv \\
        $args \\
        $mito \\
        $peak \\
        $tss \\
        $excl_regs \\
        $autosom_ref \\
        --metrics-file "${prefix}.ataqv.json" \\
        --threads $task.cpus \\
        --name $prefix \\
        $organism \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ataqv: \$( ataqv --version )
    END_VERSIONS
    """

    stub:
    def prefix        = task.ext.prefix         ?: "${meta.id}"
    def args          = task.ext.args ?: ''
    def problems_cmd  = args.contains("--log-problematic-reads") ? "echo stub > 1.problems" : ""
    """
    echo stub > ${prefix}.ataqv.json
    ${problems_cmd}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ataqv: \$( ataqv --version )
    END_VERSIONS
    """
}
