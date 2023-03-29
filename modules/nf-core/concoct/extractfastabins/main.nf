process CONCOCT_EXTRACTFASTABINS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::concoct=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/concoct:1.1.0--py38h7be5676_2':
        'quay.io/biocontainers/concoct:1.1.0--py38h7be5676_2' }"

    input:
    tuple val(meta), path(original_fasta), path(csv)

    output:
    tuple val(meta), path("${prefix}/*.fa.gz"), emit: fasta
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}

    extract_fasta_bins.py \\
        $args \\
        $original_fasta \\
        $csv \\
        --output_path ${prefix}

    ## Add prefix to each file to disambiguate one sample's 1.fa, 2.fa from sample2
    for i in ${prefix}/*.fa; do
        mv \${i} \${i/\\///${prefix}_}
        gzip \${i/\\///${prefix}_}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        concoct: \$(echo \$(concoct --version 2>&1) | sed 's/concoct //g' )
    END_VERSIONS
    """
}
