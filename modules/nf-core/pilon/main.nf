process PILON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pilon:1.24--hdfd78af_0':
        'biocontainers/pilon:1.24--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(bam), path(bai)
    val pilon_mode

    output:
    tuple val(meta), path("*.fasta") , emit: improved_assembly
    tuple val(meta), path("*.vcf")   , emit: vcf               , optional : true
    tuple val(meta), path("*.change"), emit: change_record     , optional : true
    tuple val(meta), path("*.bed")   , emit: tracks_bed        , optional : true
    tuple val(meta), path("*.wig")   , emit: tracks_wig        , optional : true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["frags", "jumps", "unpaired", "bam"]
    if ( !valid_mode.contains(pilon_mode) )  { error "Unrecognised mode to run Pilon. Options: ${valid_mode.join(', ')}" }
    def mem_mb = 1024
    if (!task.memory) {
        log.info '[Pilon] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        mem_mb = task.memory.giga < 2 ? (task.memory.mega*0.8).intValue() : task.memory.mega - 1024
    }
    """
    # `which` allows us to get the directory that contains `pilon`, independent of whether we
    # are in a container or conda environment.
    PILON_JAR=\$(dirname \$(which pilon))/../share/pilon*/pilon.jar
    java -Xmx${mem_mb}M -jar \$PILON_JAR \\
        --genome $fasta \\
        --output ${meta.id} \\
        $args \\
        --$pilon_mode $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(echo \$(pilon --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def valid_mode = ["frags", "jumps", "unpaired", "bam"]
    if ( !valid_mode.contains(pilon_mode) )  { error "Unrecognised mode to run Pilon. Options: ${valid_mode.join(', ')}" }
    """
    touch ${prefix}.fasta
    touch ${prefix}.vcf
    touch ${prefix}.change
    touch ${prefix}.bed
    touch ${prefix}.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(echo \$(pilon --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

}
