process HMMER_HMMALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1' :
        'quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1' }"

    input:
    tuple val(meta), path(fasta)
    path hmm

    output:
    tuple val(meta), path("*.sthlm.gz"), emit: sthlm
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def fastacmd = fasta.getExtension() == 'gz' ? "gunzip -c $fasta" : "cat $fasta"
    """
    $fastacmd | \\
        hmmalign \\
        $args \\
        $hmm \\
        - | gzip -c > ${meta.id}.sthlm.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmalign -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
