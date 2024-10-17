process SVDB_MERGE {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c8daa8f9d69d3c5a1a4ff08283a166c18edb0000:511069f65a53621c5503e5cfee319aa3c735abfa-0':
        'biocontainers/mulled-v2-c8daa8f9d69d3c5a1a4ff08283a166c18edb0000:511069f65a53621c5503e5cfee319aa3c735abfa-0' }"

    input:
    tuple val(meta), path(vcfs)
    val(priority)
    val(sort_inputs)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Ensure priority list matches the number of VCFs if priority is provided
    if (priority && vcfs.collect().size() != priority.collect().size()) {
        error "If priority is used, one tag per VCF is needed"
    }

    if (sort_inputs && vcfs.collect().size() > 1) {
        if (priority) {
            // make vcf-prioprity pairs and sort on VCF name, so priority is also sorted the same
            def pairs = vcfs.indices.collect { [vcfs[it], priority[it]] }
            pairs = pairs.sort { a, b -> a[0].name <=> b[0].name }
            vcfs = pairs.collect { it[0] }
            priority = pairs.collect { it[1] }
        } else {
            // if there's no priority input just sort the vcfs by name
            vcfs = vcfs.sort { it.name }
        }
    }

    // If there's only one input VCF the code above is not executed, and that VCF becomes the input
    input = vcfs

    def prio = ""
    if(priority) {
        prio = "--priority ${priority.join(',')}"
        input = ""
        for (int index = 0; index < vcfs.collect().size(); index++) {
            input += "${vcfs[index]}:${priority[index]} "
        }
    }

    """
    svdb \\
        --merge \\
        $args \\
        $prio \\
        --vcf $input \\
        > ${prefix}.vcf

    bgzip \\
        --threads ${task.cpus} \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$( echo \$(svdb) | head -1 | sed 's/usage: SVDB-\\([0-9]\\.[0-9]\\.[0-9]\\).*/\\1/' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
