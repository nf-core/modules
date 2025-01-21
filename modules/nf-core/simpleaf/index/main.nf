process SIMPLEAF_INDEX {
    tag meta.id ? "${meta.id}" : "${meta2.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.18.4--ha6fb395_1':
        'biocontainers/simpleaf:0.18.4--ha6fb395_1' }"

    input:
    tuple val(meta),  path(genome_fasta), path(genome_gtf)
    tuple val(meta2), path(transcript_fasta)

    output:
    tuple val(meta), path("${prefix}/index")                    , emit: index
    tuple val(meta), path("${prefix}/ref")                      , emit: ref, optional: true
    path "${prefix}/ref/{t2g,t2g_3col}.tsv"                     , emit: t2g, optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seq_inputs = input_args(genome_fasta, genome_gtf, transcript_fasta)//, probes_csv, features_csv)

    // Output meta needs to correspond to the input used
    meta = (transcript_fasta) ? meta2 : meta
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
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : "${meta2.id}")

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
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
    END_VERSIONS
    """
}

def input_args(genome_fasta, genome_gtf, transcript_fasta) { //, probes_csv, features_csv) {
    // if (probe_csv) {
    //     args = "--probe_csv ${probe_csv}"
    // } else if (feature_csv) {
    //     args = "--feature_csv ${feature_csv}"
    // } else
    if (transcript_fasta) {
        return "--ref-seq ${transcript_fasta}"
    } else if (genome_fasta && genome_gtf) {
        return "--fasta ${genome_fasta} --gtf ${genome_gtf}"
    } else {
        error "No valid input provided; please provide either a genome fasta + gtf set or a transcript fasta file. ${genome_fasta} ${genome_gtf} ${transcript_fasta}"
        // error "No valid input provided; please provide one of the followings: (i) a genome fasta + gtf set, (ii) a transcript fasta file, (iii) a probes csv file (iv) a features csv file."
    }

}
