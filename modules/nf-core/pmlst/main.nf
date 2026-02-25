process PMLST {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pmlst:2.0.3--hdfd78af_0':
        'biocontainers/pmlst:2.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path pmlst_db
    val scheme

    output:
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true
    tuple val(meta), path("*.txt"), emit: txt, optional: true
    tuple val("${task.process}"), val('pmlst'), eval('echo 2.0.3'), emit: versions_pmlst, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith('.gz') ? true : false
    def fasta_name = fasta.getName().replace('.gz', '')
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mkdir ${prefix}_${scheme}

    pmlst.py \\
        -i $fasta_name \\
        -o ${prefix}_${scheme}\\
        -p $pmlst_db \\
        -s ${scheme} \\
        $args

    # If output files exist, move them to expected names
    # Else touch empty files to avoid Nextflow errors
    if [ -f "${prefix}_${scheme}/results_tab.tsv" ]; then
        mv "${prefix}_${scheme}/results_tab.tsv" "${prefix}.tsv"
    else
        touch "${prefix}.tsv"
    fi

    if [ -f "${prefix}_${scheme}/results.txt" ]; then
        mv "${prefix}_${scheme}/results.txt" "${prefix}.txt"
    else
        touch "${prefix}.txt"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.txt
    """
}
