process ATLAS_SPLITMERGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::atlas=0.9.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/atlas:0.9.9--h082e891_0':
        'biocontainers/atlas:0.9.9--h082e891_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(read_group_settings), path(blacklist)

    output:
    tuple val(meta), path("*_mergedReads.bam"), path("*.txt.gz"), emit: data
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def optional = blacklist ? 'blacklist=${blacklist}' : ''
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    atlas \\
        task=splitMerge bam=${bam} \\
        readGroupSettings=${read_group_settings}\\
        $optional \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        atlas: \$((atlas 2>&1) | grep Atlas | head -n 1 | sed -e 's/^[ \t]*Atlas //')
    END_VERSIONS
    """
}
