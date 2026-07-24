process BISMARK_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bddea334e6ccbce005ce540214747acf822b040185d2198220dcfbb4b258c331/data' :
        'community.wave.seqera.io/library/bismark:3.1.0--9557d6ab108a83e4' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta, stageAs: 'tmp/*') // This change mounts as directory containing the FASTA file to prevent nested symlinks
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*bam")       , emit: bam
    tuple val(meta), path("*report.txt"), emit: report
    tuple val(meta), path("*fq.gz")     , emit: unmapped, optional: true
    tuple val("${task.process}"), val("bismark"), eval("bismark --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions, emit: versions_bismark

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if(task.ext.prefix){
        args += " --prefix ${task.ext.prefix}"
    }
    def fastq = meta.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

    // Give each strand instance `-p = cpus / n_instances` threads (one index load per strand,
    // byte-identical) instead of forking N --multicore chunks that each re-load the index.
    // minimap2/rammap can't take -p from bismark yet, so they keep the legacy --multicore path.
    def isMinimapLike = args.contains('--minimap2') || args.contains('--mm2') || args.contains('--rammap') || args.contains('--ram')
    if(isMinimapLike){

        // Legacy --multicore auto-compute (unchanged) — minimap2/rammap only.
        if(!args.contains('--multicore') && task.cpus){

            // Numbers based on recommendation by Felix for a typical mouse genome
            def ccore = 1
            def cpu_per_multicore = 3
            def mem_per_multicore = (13.GB).toBytes()
            if(args.contains('--non_directional')){
                cpu_per_multicore = 5
                mem_per_multicore = (18.GB).toBytes()
            }

            // How many multicore splits can we afford with the cpus we have?
            ccore = ((task.cpus as int) / cpu_per_multicore) as int

            // Check that we have enough memory
            try {
                def tmem = (task.memory as MemoryUnit).toBytes()
                def mcore = (tmem / mem_per_multicore) as int
                ccore = Math.min(ccore, mcore)
            } catch (Exception e) {
                log.warn "Error catched: ${e}"
                log.warn "Not able to define bismark align multicore based on available memory"
            }
            if(ccore > 1){
                args += " --multicore ${ccore}"
            }
        }
    } else if(!args.contains('--multicore') && !args.contains('--parallel') && !args.contains('-p ') && task.cpus){

        def n_instances = 2
        if(args.contains('--combined_index')){
            // combined index = one both-strands pass; parallel non-directional runs two.
            n_instances = (args.contains('--non_directional') && args.contains('--combined_index_parallel')) ? 2 : 1
        } else if(args.contains('--non_directional')){
            n_instances = 4
        }

        // bismark requires -p >= 2, so omit it below that (one thread per instance).
        def pthreads = ((task.cpus as int) / n_instances) as int
        if(pthreads >= 2){
            args += " -p ${pthreads}"
        }
    }
    """
    bismark \\
        ${fastq} \\
        --genome ${index} \\
        --bam \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.report.txt
    """
}
