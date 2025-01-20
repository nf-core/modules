process SIMPLEAF_INDEX {
    tag "$genome_fasta $transcript_fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.18.4--ha6fb395_1':
        'biocontainers/simpleaf:0.18.4--ha6fb395_1' }"

    input:
    tuple val(meta), path(genome_fasta)
    tuple val(meta2), path(genome_gtf)
    tuple val(meta3), path(transcript_fasta)

    output:
    tuple val(meta), path("${prefix}/index")                    , emit: index
    tuple val(meta), path("${prefix}/ref/{t2g,t2g_3col}.tsv")   , emit: transcript_tsv, optional: true
    tuple val(meta), path("${prefix}")                          , emit: simpleaf
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seq_inputs = (transcript_fasta) ? "--refseq $transcript_fasta" : "--gtf $genome_gtf --fasta $genome_fasta"

    // Output meta needs to correspond to the input used
    meta = (transcript_fasta) ? meta3 : meta
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # export required var
    export ALEVIN_FRY_HOME=.

    # set maximum number of file descriptors for temp files
    ulimit -n 2048

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf index
    simpleaf \\
        index \\
        --threads $task.cpus \\
        $seq_inputs \\
        $args \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevin-fry: \$(alevin-fry --version | sed -e "s/alevin-fry //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : "${meta3.id}")

    """
    mkdir -p ${prefix}/index
    mkdir -p ${prefix}/ref
    touch ${prefix}/index/piscem_idx_cfish.json
    touch ${prefix}/index/piscem_idx.ectab
    touch ${prefix}/index/piscem_idx.sshash
    touch ${prefix}/ref/t2g_3col.tsv
    touch ${prefix}/ref/roers_ref.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevin-fry: \$(alevin-fry --version | sed -e "s/alevin-fry //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
    END_VERSIONS
    """
}
