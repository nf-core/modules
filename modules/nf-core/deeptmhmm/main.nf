process DEEPTMHMM {
    tag '$deeptmhmm'
    label 'process_medium'

    conda "bioconda::pybiolib=1.1.1367"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pybiolib:1.1.1367--pyhdfd78af_0':
        'biocontainers/pybiolib:1.1.1367--pyhdfd78af_0' }"

    input:
    path fasta

    output:
    path "biolib_results/TMRs.gff3"                 , emit: gff3
    path "biolib_results/predicted_topologies.3line", emit: line3
    path "biolib_results/deeptmhmm_results.md"      , emit: md
    path "biolib_results/*_probs.csv"               , optional: true, emit: csv
    path "biolib_results/plot.png"                  , optional: true, emit: png
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed = fasta.name.endsWith(".gz")
    def fasta_name = fasta.name.replace(".gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    biolib \\
        run \\
        DTU/DeepTMHMM \\
        --fasta ${fasta_name} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biolib: \$(echo \$(biolib --version) | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    mkdir biolib_results
    touch biolib_results/TMRs.gff3
    touch biolib_results/predicted_topologies.3line
    touch biolib_results/deeptmhmm_results.md
    touch biolib_results/MX_probs.csv
    touch biolib_results/plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biolib: \$(echo \$(biolib --version) | sed -n 's/.*version \\([0-9.]*\\).*/\\1/p' )
    END_VERSIONS
    """
}
