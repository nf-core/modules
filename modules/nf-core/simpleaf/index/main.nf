process SIMPLEAF_INDEX {
    tag "$genome_fasta $transcript_fasta"
    label 'process_high'

    conda "bioconda::simpleaf=0.14.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.14.1--h4ac6f70_0':
        'biocontainers/simpleaf:0.14.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(genome_fasta)
    tuple val(meta2), path(genome_gtf)
    path transcript_fasta

    output:
    tuple val(meta), path("salmon/index")              , emit: index
    tuple val(meta), path("salmon/ref/t2g_3col.tsv")   , emit: transcript_tsv, optional: true
    tuple val(meta), path("salmon")                    , emit: salmon
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seq_inputs = (transcript_fasta) ? "--refseq $transcript_fasta" : "--gtf $genome_gtf --fasta $genome_fasta"
    """
    # export required var
    export ALEVIN_FRY_HOME=.

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
        simpleaf: \$(simpleaf -V | tr -d '\\n' | cut -d ' ' -f 2)
        salmon: \$(salmon --version | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
