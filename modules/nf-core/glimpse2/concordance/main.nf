process GLIMPSE2_CONCORDANCE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::glimpse-bio=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.0--hf340a29_0':
        'quay.io/biocontainers/glimpse-bio:2.0.0--hf340a29_0' }"

    input:
    tuple val(meta), path(estimate), path(estimate_index), path(truth), path(truth_index), path(freq), path(freq_index), path(samples), val(region)
    tuple val(meta2), path(groups), val(bins), val(ac_bins), val(allele_counts)
    val(min_val_gl)
    val(min_val_dp)

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
    def samples_cmd  = samples         ? "--samples ${samples}"             : ""
    def groups_cmd   = groups          ? "--groups ${groups}"               : ""
    def bins_cmd     = bins            ? "--bins ${bins}"                   : ""
    def ac_bins_cmd  = ac_bins         ? "--ac-bins ${ac_bins}"             : ""
    def ale_ct_cmd   = allele_counts   ? "--allele-counts ${allele_counts}" : ""

    if (((groups ? 1:0) + (bins ? 1:0) + (ac_bins ? 1:0) + (allele_counts ? 1:0)) != 1) error "One and only one argument should be selected between groups, bins, ac_bins, allele_counts"

    """
    echo $region $freq $truth $estimate > input.txt
    GLIMPSE2_concordance \\
        $args \\
        $samples_cmd \\
        $groups_cmd \\
        $bins_cmd \\
        $ac_bins_cmd \\
        $ale_ct_cmd \\
        --min-val-gl $min_val_gl \\
        --min-val-dp $min_val_dp \\
        --input input.txt \\
        --thread $task.cpus \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            glimpse2: "\$(GLIMPSE2_concordance --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
