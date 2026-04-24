process DIAMOND_BLASTX {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/diamond:2.1.24--hf93d47f_0'
        : 'quay.io/biocontainers/diamond:2.1.24--hf93d47f_0'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(db)
    val out_ext
    val blast_columns

    output:
    tuple val(meta), path('*.blast'), optional: true, emit: blast
    tuple val(meta), path('*.xml'), optional: true, emit: xml
    tuple val(meta), path('*.txt'), optional: true, emit: txt
    tuple val(meta), path('*.daa'), optional: true, emit: daa
    tuple val(meta), path('*.sam'), optional: true, emit: sam
    tuple val(meta), path('*.tsv'), optional: true, emit: tsv
    tuple val(meta), path('*.paf'), optional: true, emit: paf
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('diamond'), eval("diamond --version | sed 's/diamond version //g'"), emit: versions_diamond, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    def columns = blast_columns ? "${blast_columns}" : ''
    if (out_ext == 'blast') {
        outfmt = 0
    }
    else if (out_ext == 'xml') {
        outfmt = 5
    }
    else if (out_ext == 'txt') {
        outfmt = 6
    }
    else if (out_ext == 'daa') {
        outfmt = 100
    }
    else if (out_ext == 'sam') {
        outfmt = 101
    }
    else if (out_ext == 'tsv') {
        outfmt = 102
    }
    else if (out_ext == 'paf') {
        outfmt = 103
    }
    else {
        outfmt = 6
        out_ext = 'txt'
        log.warn("Unknown output file format provided (${out_ext}): selecting DIAMOND default of tabular BLAST output (txt)")
    }
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    DB=`find -L ./ -name "*.dmnd" | sed 's/\\.dmnd\$//'`

    diamond \\
        blastx \\
        --threads ${task.cpus} \\
        --db \$DB \\
        --query ${fasta_name} \\
        --outfmt ${outfmt} ${columns} \\
        ${args} \\
        --out ${prefix}.${out_ext} \\
        --log

    mv diamond.log ${prefix}.log
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (out_ext == 'blast') {
        outfmt = 0
    }
    else if (out_ext == 'xml') {
        outfmt = 5
    }
    else if (out_ext == 'txt') {
        outfmt = 6
    }
    else if (out_ext == 'daa') {
        outfmt = 100
    }
    else if (out_ext == 'sam') {
        outfmt = 101
    }
    else if (out_ext == 'tsv') {
        outfmt = 102
    }
    else if (out_ext == 'paf') {
        outfmt = 103
    }
    else {
        outfmt = 6
        out_ext = 'txt'
        log.warn("Unknown output file format provided (${out_ext}): selecting DIAMOND default of tabular BLAST output (txt)")
    }

    """
    echo "${args}"
    touch ${prefix}.${out_ext}
    touch ${prefix}.log
    """
}
