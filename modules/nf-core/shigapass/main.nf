process SHIGAPASS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/shigapass:1.5.0--hdfd78af_0':
        'biocontainers/shigapass:1.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.tsv"),                 emit: report
    tuple val(meta), path("*_ShigaPass_Flex_summary.tsv"),  emit: flex_tsv, optional: true
    tuple val("${task.process}"), val('shigapass'), eval("ShigaPass.sh -v 2>&1 | sed 's/^.*ShigaPass version //'"), topic: versions, emit: versions_shigapass

    when:
        task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def dbPath = (
        workflow.containerEngine == 'singularity' ||
        workflow.containerEngine == 'apptainer'
    ) ?
        "/usr/local/share/shigapass-1.5.0/db" :
        "\$CONDA_PREFIX/share/shigapass-1.5.0/db"
    def is_compressed = fasta.getName().endsWith(".gz")
    def fasta_name    = fasta.getName().replace(".gz", "")
    """
    # Optional decompression
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # Create temp file
    ls ${fasta_name} > ${fasta_name}_tmp.txt

    ShigaPass.sh \\
        -l ${fasta_name}_tmp.txt \\
        $args \\
        -p ${dbPath} \\
        -t ${task.cpus} \\
        -o ${prefix}

    # Remove temp file
    rm ${fasta_name}_tmp.txt

    # Change from ; delim to tab delim
    sed 's/;/\\t/g' ${prefix}/ShigaPass_summary.csv > ${prefix}.tsv
    [ ! -f ${prefix}/ShigaPass_Flex_summary.csv ] || sed 's/;/\\t/g' ${prefix}/ShigaPass_Flex_summary.csv > ${prefix}_Flex_summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
