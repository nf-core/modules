process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.21--h5bf99c6_0' :
        'quay.io/biocontainers/minimap2:2.21--h5bf99c6_0' }"

    input:
    tuple val(meta), path(reads)
    path reference
    val sam_format
    val preset_pacbio_reads
    val preset_nanopore_reads
    val preset_pacbio_hifi_reads
    val preset_pacbio_overlap
    val preset_nanopore_overlap
    val preset_asm5
    val preset_asm10
    val preset_asm20
    val preset_nanopore_spliced
    val preset_pacbio_spliced
    val preset_short_read

    output:
    tuple val(meta), path("*.paf"), emit: paf, optional: true
    tuple val(meta), path("*.sam"), emit: sam, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "$reads" : "${reads[0]} ${reads[1]}"
    def sam_output = sam_format ? "-a -o ${prefix}.sam" : "-o ${prefix}.paf"
    def preset = preset_pacbio_reads ? "-x map-pb" : preset_nanopore_reads ? "-x map-ont" : preset_pacbio_hifi_reads ? "-x map-hifi" :\
                 preset_pacbio_overlap ? "-x ava-pb" : preset_nanopore_overlap ? "-x ava-ont" : preset_asm5 ? "-x asm5" : preset_asm10 ? "-x asm10" :\
                 preset_asm20 ? "-x asm20" : preset_nanopore_spliced ? "-x splice" : preset_pacbio_spliced ? "-x splice:hq" : preset_short_read ? "-x sr" : ''
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $reference \\
        $input_reads \\
        $sam_output \\
        $preset


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
