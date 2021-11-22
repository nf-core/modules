process HMMER_HMMALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1"
    } else {
        container "quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1"
    }

    input:
    tuple val(meta), path(fasta)
    path hmm

    output:
    tuple val(meta), path("*.sthlm.gz"), emit: sthlm
    path "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fastacmd = fasta.getExtension() == 'gz' ? "gunzip -c $fasta" : "cat $fasta"
    """
    $fastacmd | \\
        hmmalign \\
        $options.args \\
        $hmm \\
        - | gzip -c > ${meta.id}.sthlm.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(hmmalign -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}
