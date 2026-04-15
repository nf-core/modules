process CONCOCT_EXTRACTFASTABINS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py39h8907335_8'
        : 'biocontainers/concoct:1.1.0--py39h8907335_8'}"

    input:
    tuple val(meta), path(original_fasta), path(csv)

    output:
    tuple val(meta), path("${prefix}/*.fa.gz"), emit: fasta
    tuple val("${task.process}"), val('concoct'), eval("concoct --version 2>&1 | sed -n 's/concoct //p'"), topic: versions, emit: versions_concoct

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    extract_fasta_bins.py \\
        ${args} \\
        ${original_fasta} \\
        ${csv} \\
        --output_path ${prefix}

    ## Add prefix to each file to disambiguate one sample's 1.fa, 2.fa from sample2
    for i in ${prefix}/*.fa; do
        mv \${i} \${i/\\///${prefix}_}
        gzip \${i/\\///${prefix}_}
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    echo "" | gzip > ${prefix}/${prefix}.fa.gz
    """
}
