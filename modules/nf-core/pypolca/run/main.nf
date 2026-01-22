process PYPOLCA_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pypolca:0.4.0--pyhdfd78af_0':
        'biocontainers/pypolca:0.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(contigs)

    output:
    tuple val(meta), path("${prefix}/*_corrected.fasta"), emit: polished
    tuple val(meta), path("${prefix}/*.vcf")            , emit: vcf
    tuple val(meta), path("${prefix}/*.report")         , emit: report
    tuple val("${task.process}"), val('pypolca'), eval('pypolca --version'), emit: versions_pypolca, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: '--careful'
    prefix            = task.ext.prefix ?: "${meta.id}"
    def read_files    = reads instanceof List ? reads : [reads]
    def read_file_arg = read_files.size() > 1 ? "-1 ${read_files[0]} -2 ${read_files[1]}" : "-1 ${read_files[0]}"
    """
    gzip -cdf $contigs > contigs_uncompressed

    pypolca \
        run \
        -a contigs_uncompressed \\
        $read_file_arg \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        --prefix ${prefix} \\
        $args
    """

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    mkdir $prefix
    touch $prefix/${prefix}_corrected.fasta
    touch $prefix/${prefix}.vcf
    touch $prefix/${prefix}.report
    """
}
