process SVDB_MERGE {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c8daa8f9d69d3c5a1a4ff08283a166c18edb0000:511069f65a53621c5503e5cfee319aa3c735abfa-0':
        'biocontainers/mulled-v2-c8daa8f9d69d3c5a1a4ff08283a166c18edb0000:511069f65a53621c5503e5cfee319aa3c735abfa-0' }"

    input:
    tuple val(meta), path(vcfs)
    val(input_priority)
    val(sort_inputs)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    println vcfs.collect()
    println input_priority.collect().size()

    // Ensure priority list matches the number of VCFs if priority is provided
    if (input_priority && vcfs.collect().size() != input_priority.collect().size()) {
        error "If priority is used, one tag per VCF is needed"
    }

    def input = ""
    def prio = ""
    if (input_priority) {
        if (vcfs.collect().size() > 1 && sort_inputs) {
            // make vcf-prioprity pairs and sort on VCF name, so priority is also sorted the same
            def pairs = vcfs.indices.collect { [vcfs[it], input_priority[it]] }
            pairs = pairs.sort { a, b -> a[0].name <=> b[0].name }
            vcfs = pairs.collect { it[0] }
            priority = pairs.collect { it[1] }
        } else {
            priority = input_priority
        }

        // Build inputs
        prio = "--priority ${input_priority.join(',')}"
        input = ""
        for (int index = 0; index < vcfs.collect().size(); index++) {
            input += "${vcfs[index]}:${priority[index]} "
        }

    } else {
        // if there's no priority input just sort the vcfs by name if possible
        input = (vcfs.collect().size() > 1 && sort_inputs) ? vcfs.sort { it.name } : vcfs
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
