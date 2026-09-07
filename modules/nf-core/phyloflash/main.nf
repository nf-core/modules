process PHYLOFLASH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d993344c3f636cb0cca9519b11fcf30faafb13fca2fd33090104e5f52d8fd643/data' :
        'community.wave.seqera.io/library/phyloflash:3.4.2--87628969a9477d43' }"

    input:
    tuple val(meta), path(reads)
    path  silva_db
    path  univec_db

    output:
    tuple val(meta), path("${meta.id}*/*"), emit: results
    tuple val("${task.process}"), val('phyloflash'), eval("phyloFlash.pl -version 2>&1 | sed 's/^.*phyloFlash v//'"), topic: versions, emit: versions_phyloflash

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        phyloFlash.pl \\
            ${args} \\
            -read1 ${reads[0]} \\
            -lib ${prefix} \\
            -interleaved \\
            -dbhome . \\
            -CPUs ${task.cpus}

        mkdir ${prefix}
        mv ${prefix}.* ${prefix}
        """
    } else {
        """
        phyloFlash.pl \\
            ${args} \\
            -read1 ${reads[0]} \\
            -read2 ${reads[1]} \\
            -lib ${prefix} \\
            -dbhome . \\
            -CPUs ${task.cpus}

        mkdir ${prefix}
        mv ${prefix}.* ${prefix}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.SSU.collection.fasta
    touch ${prefix}/${prefix}.phyloFlash
    """
}
