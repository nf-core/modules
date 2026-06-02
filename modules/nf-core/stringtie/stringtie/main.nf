process STRINGTIE_STRINGTIE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3facd74a0f728c9bb9e9a731b58c343895d2dbdfeb812ce5747f701103fc61cf/data' :
        'community.wave.seqera.io/library/stringtie:3.0.3--e8043d00caecd051' }"

    input:
    tuple val(meta), path(srbam), path(lrbam)
    path(annotation_gtf)
    val(mode)

    output:
    tuple val(meta), path("${prefix}.ballgown/")         , optional: true, emit: ballgown
    tuple val(meta), path("${prefix}.coverage.gtf")      , optional: true, emit: coverage_gtf
    tuple val(meta), path("${prefix}.transcripts.gtf")   , emit: transcript_gtf
    tuple val(meta), path("${prefix}.gene.abundance.txt"), emit: abundance
    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version'), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def args2     = task.ext.args2 ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"

    def reference = annotation_gtf ? "-G $annotation_gtf" : ""
    def ballgown  = annotation_gtf ? "-b ${prefix}.ballgown" : ""
    def coverage  = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""

    def run_mode   = ''
    if ("${mode}") {
        def valid_run_modes = ['expression-estimation', 'long-reads-assembly', 'mix-reads-assembly']
        if (!(mode in valid_run_modes)) {
            error "Invalid mode: ${mode}. Valid options are: ${valid_run_modes.join(', ')}"
        }
        run_mode  = "${mode}" ? (
            mode == 'expression-estimation' ? '-L -e':
            mode == 'long-reads-assembly' ? '-L' :
            mode == 'mix-reads-assembly' ? '--mix' : ''
        ): ''
    }

    def bam_inputs = ''
    if (srbam && !lrbam) {
        bam_inputs = "$srbam "
    } else if (!srbam && lrbam) {
        bam_inputs = "$lrbam "
    } else if (srbam && lrbam) {
        bam_inputs = (mode == 'mix-reads-assembly') ? "$srbam $lrbam" : "$srbam $lrbam"
    }
    
    """
    stringtie \\
        -o ${prefix}.transcripts.gtf \\
        ${run_mode} \\
        ${args2} \\
        ${reference} \\
        -A ${prefix}.gene.abundance.txt \\
        ${coverage} \\
        ${ballgown} \\
        -p ${task.cpus} \\
        ${bam_inputs} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}.ballgown/
    touch ${prefix}.transcripts.gtf
    touch ${prefix}.gene.abundance.txt
    touch ${prefix}.coverage.gtf
    """
}