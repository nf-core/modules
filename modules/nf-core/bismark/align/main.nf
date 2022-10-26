process BISMARK_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bismark=0.23.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.23.0--0' :
        'quay.io/biocontainers/bismark:0.23.0--0' }"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*bam")       , emit: bam
    tuple val(meta), path("*report.txt"), emit: report
    tuple val(meta), path("*fq.gz")     , optional:true, emit: unmapped
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Try to assign sensible bismark memory units according to what the task was given
    def ccore = 1
    def cpu_per_multicore = 3
    def mem_per_multicore = (13.GB).toBytes()

    if( task.cpus ){
        // Numbers based on recommendation by Felix for a typical mouse genome
        if( params.single_cell || params.zymo || params.non_directional ){
            cpu_per_multicore = 5
            mem_per_multicore = (18.GB).toBytes()
        }
        // Check if the user has specified this and overwrite if so
        if(params.bismark_align_cpu_per_multicore) {
            cpu_per_multicore = (params.bismark_align_cpu_per_multicore as int)
        }
        if(params.bismark_align_mem_per_multicore) {
            mem_per_multicore = (params.bismark_align_mem_per_multicore as nextflow.util.MemoryUnit).toBytes()
        }
        // How many multicore splits can we afford with the cpus we have?
        ccore = ((task.cpus as int) / cpu_per_multicore) as int
        // Check that we have enough memory, assuming 13GB memory per instance (typical for mouse alignment)
        try {
            def tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
            def mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (all) {
            log.warn "Not able to define bismark align multicore based on available memory"
        }
    }
    def multicore  = (ccore > 1) ? "--multicore ${ccore}" : ''

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq      = meta.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    bismark \\
        $fastq \\
        $args \\
        $multicore \\
        --genome $index \\
        --bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
