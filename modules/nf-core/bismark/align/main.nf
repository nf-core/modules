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

    // Aligner thread model:
    //  - Bowtie 2 + HISAT2 run 2 (directional/pbat) or 4 (non-directional) *concurrent*
    //    strand instances, each of which honours `-p` (threads per instance). Dividing the
    //    task's cpus across the instances (`-p = task.cpus / n_instances`) uses all cpus with
    //    a single index load per instance, instead of `--multicore N` forking N chunks that
    //    each re-load the index (Bowtie 2: N x 2/4 concurrent index copies -> high peak RAM).
    //    Bowtie 2 `-p N --reorder` is thread-invariant (byte-identical); HISAT2 is near-
    //    invariant (a handful of reads re-baseline). See the resource-model benchmark.
    //  - minimap2/rammap ignore `-p` (bismark hard-codes minimap2 `-t 2`; rammap auto-scales),
    //    so they keep the legacy `--multicore` auto-compute untouched.
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

        // Bowtie 2 / HISAT2: split the cpus across the concurrent strand instances.
        def n_instances = 2
        if(args.contains('--combined_index')){
            // combined = one both-strands pass (n=1); only the opt-in parallel non-directional
            // path runs two concurrent passes (n=2).
            n_instances = (args.contains('--non_directional') && args.contains('--combined_index_parallel')) ? 2 : 1
        } else if(args.contains('--non_directional')){
            n_instances = 4
        }

        // -p threads per instance; emit only when it lands >= 2 (bismark requires -p >= 2),
        // otherwise omit it (faithful single-thread-per-instance default).
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
