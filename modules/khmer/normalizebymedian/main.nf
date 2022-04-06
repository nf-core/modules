process KHMER_NORMALIZEBYMEDIAN {
    tag "${name}"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    path pe_reads
    path se_reads
    val  name

    output:
    path "${name}.fastq.gz", emit: reads
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    pe_args = pe_reads ? "--paired" : ""
    se_args = se_reads ? "--unpaired-reads ${se_reads}" : ""
    files   = pe_reads ? pe_reads : se_reads
    """
    normalize-by-median.py \\
        -M ${task.memory.toGiga()}e9 \\
        --gzip $args \\
        -o ${name}.fastq.gz \\
        $pe_args \\
        $se_args \\
        $files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
