process GSTAMA_POLYACLEANUP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0':
        'quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_tama.fa.gz")                   , emit: fasta
    tuple val(meta), path("*_tama_polya_flnc_report.txt.gz"), emit: report
    tuple val(meta), path("*_tama_tails.fa.gz")             , emit: tails
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$fasta" == "${prefix}.fasta" | "$fasta" == "${prefix}.fa" ) error "Input and output names are the same, set prefix in module configuration"
    """
    tama_flnc_polya_cleanup.py \\
        -f $fasta \\
        -p ${prefix} \\
        $args
    gzip ${prefix}.fa
    gzip ${prefix}_polya_flnc_report.txt
    gzip ${prefix}_tails.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gstama: \$( tama_collapse.py -version | grep 'tc_version_date_'|sed 's/tc_version_date_//g' )
    END_VERSIONS
    """
}
