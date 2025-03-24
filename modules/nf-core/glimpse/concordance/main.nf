process GLIMPSE_CONCORDANCE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:1.1.1--hce55b13_1':
        'biocontainers/glimpse-bio:1.1.1--hce55b13_1' }"

    input:
    tuple val(meta), path(estimate), path(estimate_index), path(freq), path(freq_index), path(truth), path(truth_index), val(region)
    val(min_prob)
    val(min_dp)
    val(bins)

    output:
    tuple val(meta), path("*.error.cal.txt.gz")  , emit: errors_cal
    tuple val(meta), path("*.error.grp.txt.gz")  , emit: errors_grp
    tuple val(meta), path("*.error.spl.txt.gz")  , emit: errors_spl
    tuple val(meta), path("*.rsquare.grp.txt.gz"), emit: rsquare_grp
    tuple val(meta), path("*.rsquare.spl.txt.gz"), emit: rsquare_spl
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def min_prob_cmd = min_prob ? "--minPROB ${min_prob}" : "--minPROB 0.9999"
    def min_dp_cmd   = min_dp   ? "--minDP ${min_dp}"     : "--minDP 8"
    def bins_cmd     = bins     ? "--bins ${bins}"        : "--bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000"
    """
    echo $region $freq $truth $estimate > input.txt
    GLIMPSE_concordance \\
        $args \\
        --input input.txt \\
        --thread $task.cpus \\
        --output ${prefix} \\
        $min_prob_cmd \\
        $min_dp_cmd \\
        $bins_cmd

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            glimpse: "\$(GLIMPSE_concordance --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def args    = task.ext.args   ?: ""
    """
    echo "" | gzip > ${prefix}.error.cal.txt.gz
    echo "" | gzip > ${prefix}.error.grp.txt.gz
    echo "" | gzip > ${prefix}.error.spl.txt.gz
    echo "" | gzip > ${prefix}.rsquare.grp.txt.gz
    echo "" | gzip > ${prefix}.rsquare.spl.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glimpse: "\$(GLIMPSE_concordance --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]')"
    END_VERSIONS
    """
}
