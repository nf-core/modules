process RIBOCODE_PREPARE {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribocode:1.2.15--pyhfa5458b_0':
        'biocontainers/ribocode:1.2.15--pyhfa5458b_0' }"

    input:
    path fasta
    path gtf

    output:
    path "annotation"             , emit: annotation
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # nf-core: ensure FASTA is uncompressed
    GENOME="$fasta"
    if [[ "\$GENOME" == *.gz ]]; then
        gunzip -c "\$GENOME" > genome.fa
        GENOME="genome.fa"
    fi

    # Update the GTF
    GTFupdate \\
        $gtf  \\
        > updated.gtf

    # Run prepare_transcripts with uncompressed FASTA
    prepare_transcripts \\
        -g updated.gtf \\
        -f "\$GENOME" \\
        -o annotation

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RiboCode: \$(RiboCode --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir annotation

    touch annotation/transcripts_cds.txt
    touch annotation/transcripts_sequence.fa
    touch annotation/transcripts.pickle

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RiboCode: \$(echo "1.2.15")
    END_VERSIONS
    """
}
