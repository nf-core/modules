process GSTAMA_POLYACLEANUP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gs-tama:1.0.3--hdfd78af_0':
        'quay.io/biocontainers/gs-tama:1.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_seq.fa.gz")               , emit: seq
    tuple val(meta), path("*_polya_flnc_report.txt.gz"), emit: report
    tuple val(meta), path("*_tails.fa.gz")             , emit: tails
    tuple val("${task.process}"), val('gstama'), eval("tama_collapse.py -version | sed -n 's/tc_version_date_//p'"), emit: versions_gstama, topic: versions

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

    mv ${prefix}.fa ${prefix}_seq.fa

    gzip ${prefix}_seq.fa
    gzip ${prefix}_polya_flnc_report.txt
    gzip ${prefix}_tails.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_seq.fa.gz
    echo "" | gzip > ${prefix}_polya_flnc_report.txt.gz
    echo "" | gzip > ${prefix}_tails.fa.gz
    """
}
