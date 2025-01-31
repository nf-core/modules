process LEARNMSA_ALIGN {
    tag "$meta.id"
    label 'process_medium'
    container "registry.hub.docker.com/felbecker/learnmsa:2.0.9"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln{.gz,}"), emit: alignment
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("LearnMSA align module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    """
    learnMSA \\
        -i $fasta \\
        -o "${prefix}.aln" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        learnmsa: \$(learnMSA -h | grep 'version' | awk -F 'version ' '{print \$2}' | awk '{print \$1}' | sed 's/)//g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        learnmsa: \$(if command -v learnMSA &>/dev/null; then learnMSA -h | grep 'version' | awk -F 'version ' '{print \$2}' | awk '{print \$1}' | sed 's/)//g'; else echo "STUB_TEST_HARDCODED_VERSION"; fi)
    END_VERSIONS
    """
}
