process BISMARK_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bismark=0.24.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.0--hdfd78af_0' :
        'biocontainers/bismark:0.24.0--hdfd78af_0' }"

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
    def args = task.ext.args ?: ''
    if(task.ext.prefix){
        args += " --prefix ${task.ext.prefix}"
    }
    def fastq = meta.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

    // Try to assign sensible bismark --multicore if not already set
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
            def tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
            def mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (all) {
            log.warn "Not able to define bismark align multicore based on available memory"
        }
        if(ccore > 1){
            args += " --multicore ${ccore}"
        }
    }

    """
    bismark \\
        $fastq \\
        --genome $index \\
        --bam \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
