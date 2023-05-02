process NEXTDENOVO {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::nextdenovo=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nextdenovo%3A2.5.2--py39h20169af_0':
        'quay.io/biocontainers/nextdenovo:2.5.2--py39h20169af_0' }"

    input:
    tuple val(meta), path(reads)
    val (read_type) // clr, ont, hifi
    val (genome_size) // estimated genome size, suffix K/M/G recognized

    output:
    tuple val(meta), path("*cns.fasta.gz")  , emit: corrected_reads
    tuple val(meta), path("*asm.fasta.gz")  , emit: assembly
    tuple val(meta), path("*asm.fasta.stat"), emit: stat
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #Prepare the input file
    ls ${reads} > input.fofn

    #Prepare the conf file run.cfg
    echo "[General]
        job_type = local
        job_prefix = ${prefix}
        task = all
        rewrite = yes
        deltmp = yes
        parallel_jobs = 4
        input_type = raw
        read_type = ${read_type}
        input_fofn = input.fofn
        workdir = $prefix

        [correct_option]
        read_cutoff = 1k
        genome_size = $genome_size
        sort_options = -t ${task.cpus}
        minimap2_options_raw = -t ${task.cpus}
        pa_correction = 1
        correction_options = -p 10

        [assemble_option]
        minimap2_options_cns = -t ${task.cpus}
        nextgraph_options = -a 1" > run.cfg

    #Run
    nextDenovo run.cfg

    #Data organization
    gzip -c ${prefix}/02.cns_align/01.seed_cns.sh.work/seed_cns*/cns.fasta > ${prefix}_cns.fasta.gz
    gzip -c ${prefix}/03.ctg_graph/nd.asm.fasta > ${prefix}_nd.asm.fasta.gz
    mv ${prefix}/03.ctg_graph/nd.asm.fasta.stat ${prefix}_nd.asm.fasta.stat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextDenovo : \$(echo \$(nextDenovo --version 2>&1) | sed 's/^.*nextDenovo //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
