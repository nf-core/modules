process GATK4_ANNOTATEINTERVALS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(mappable_regions)
    tuple val(meta6), path(mappable_regions_tbi)
    tuple val(meta7), path(segmental_duplication_regions)
    tuple val(meta8), path(segmental_duplication_regions_tbi)

    output:
    tuple val(meta), path("*.tsv"), emit: annotated_intervals
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def inputs = intervals.collect(){ "--intervals ${it}" }.join(" ")
    def mappability_track = mappable_regions ? "--mappability-track ${mappable_regions}" : ""
    def segmental_duplication_tracks = segmental_duplication_regions ? "--segmental-duplication-track ${segmental_duplication_regions}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK AnnotateIntervals] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        AnnotateIntervals \\
        ${inputs} \\
        --reference ${fasta} \\
        --output ${prefix}.tsv \\
        ${mappability_track} \\
        ${segmental_duplication_tracks} \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
