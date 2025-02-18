// NOTE because the default indexer, piscem, needs to frequently read and write a large number of intermediate files, if your use case involves the situations where the CPU and storage are not physically connected, we recommend setting `--work-dir /path/to/a/local/dir` or in the `ext.args` in nextflow.config, or  `scratch = true`, to avoid runtime issues.
process SIMPLEAF_INDEX {
    tag "${meta.id ?: meta2.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.19.2--ha6fb395_0':
        'biocontainers/simpleaf:0.19.2--ha6fb395_0' }"

    input:
    tuple val(meta),  path(genome_fasta), path(genome_gtf)
    tuple val(meta2), path(transcript_fasta)
    tuple val(meta3), path(probe_csv)
    tuple val(meta4), path(feature_csv)

    output:
    tuple val(meta), path("${prefix}/index")                    , emit: index
    tuple val(meta), path("${prefix}/ref")                      , emit: ref, optional: true
    tuple val(meta), path("${prefix}/ref/{t2g,t2g_3col}.tsv")   , emit: t2g, optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    (meta, seq_inputs) = input_args(genome_fasta, genome_gtf, transcript_fasta, probe_csv, feature_csv, meta, meta2, meta3, meta4)

    // Output meta needs to correspond to the input used
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
    meta = meta ? meta : [id: 'stub']
    prefix = task.ext.prefix ?: "${meta.id}"

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

def input_args(genome_fasta, genome_gtf, transcript_fasta, probe_csv, feature_csv, meta, meta2, meta3, meta4) {
    // check if all null
    if (!genome_fasta && !genome_gtf && !transcript_fasta && !probe_csv && !feature_csv) {
        error "No valid input provided; please provide either a genome fasta + gtf set or a transcript fasta file."
    }

    if (feature_csv) {
        return [meta4, "--feature-csv ${feature_csv}"]
    } else if (probe_csv) {
        return [meta3, "--probe-csv ${probe_csv}"]
    } else if (transcript_fasta) {
        return [meta2, "--ref-seq ${transcript_fasta}"]
    } else if (genome_fasta && genome_gtf) {
        return [meta, "--fasta ${genome_fasta} --gtf ${genome_gtf}"]
    } else {
        error "No valid input provided; please provide one of the followings: (i) a genome fasta + gtf set, (ii) a transcript fasta file, (iii) a probes csv file (iv) a features csv file."
    }

}
