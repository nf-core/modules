process GLIMPSE2_CONCORDANCE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::glimpse-bio=2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/glimpse-bio:2.0.0--hf340a29_0':
        'quay.io/biocontainers/glimpse-bio:2.0.0--hf340a29_0' }"

    input:
    tuple val(meta), val(region), path(freq), path(truth), path(estimate), path(samples)
    tuple val(meta2), path(groups)

    output:
    tuple val(meta), path("*.")  , emit: errors
    tuple val(meta), path("*.rsquare.*.txt.gz"), emit: rsquare
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def samples_cmd  = samples         ? "--samples ${samples}" : ""
    def groups_cmd   = groups          ? "--groups ${groups}"   : ""
    """
    echo $region $freq $truth $estimate > input.txt
    GLIMPSE2_concordance \\
        $args \\
        $samples_cmd \\
        $groups_cmd \\
        --input input.txt \\
        --thread $task.cpus \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            glimpse2: "\$(GLIMPSE2_concordance --help | sed -nr '/Version/p' | grep -o -E '([0-9]+.){1,2}[0-9]' | head -1)"
    END_VERSIONS
    """
}
