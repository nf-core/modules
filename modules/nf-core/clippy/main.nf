process CLIPPY {
    tag "$meta.id"
    label "process_high"

    conda "bioconda::clippy=1.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clippy:1.5.0--pyhdfd78af_0' :
        'biocontainers/clippy:1.5.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bed)
    path gtf
    path fai

    output:
    tuple val(meta), path("*_Peaks.bed")             ,emit: peaks
    tuple val(meta), path("*_Summits.bed")           ,emit: summits
    tuple val(meta), path("*_intergenic_regions.gtf"),emit: intergenic_gtf, optional: true
    path "versions.yml"                              ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    clippy -i $bed \
        -o $prefix \
        -a $gtf \
        -g $fai \
        -t ${task.cpus} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clippy: \$(clippy -v)
    END_VERSIONS
    """
}
