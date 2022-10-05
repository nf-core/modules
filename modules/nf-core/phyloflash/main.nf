process PHYLOFLASH {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::phyloflash=3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phyloflash:3.4--hdfd78af_1' :
        'quay.io/biocontainers/phyloflash:3.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    path  silva_db
    path  univec_db

    output:
    tuple val(meta), path("${meta.id}*/*"), emit: results
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        phyloFlash.pl \\
            $args \\
            -read1 ${reads[0]} \\
            -lib $prefix \\
            -interleaved \\
            -dbhome . \\
            -CPUs $task.cpus

        mkdir $prefix
        mv ${prefix}.* $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            phyloflash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
        END_VERSIONS
        """
    } else {
        """
        phyloFlash.pl \\
            $args \\
            -read1 ${reads[0]} \\
            -read2 ${reads[1]} \\
            -lib $prefix \\
            -dbhome . \\
            -CPUs $task.cpus

        mkdir $prefix
        mv ${prefix}.* $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            phyloflash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.SSU.collection.fasta
    touch ${prefix}/${prefix}.phyloFlash

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyloflash: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
