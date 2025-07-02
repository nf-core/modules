process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/adapterremoval:2.3.2--hb7ba0dd_0' :
        'biocontainers/adapterremoval:2.3.2--hb7ba0dd_0' }"

    input:
    tuple val(meta), path(reads)
    path(adapterlist)

    output:
    tuple val(meta), path("${prefix}.truncated.fastq.gz")          , emit: singles_truncated  , optional: true
    tuple val(meta), path("${prefix}.discarded.fastq.gz")          , emit: discarded          , optional: true
    tuple val(meta), path("${prefix}.pair{1,2}.truncated.fastq.gz"), emit: paired_truncated   , optional: true
    tuple val(meta), path("${prefix}.collapsed.fastq.gz")          , emit: collapsed          , optional: true
    tuple val(meta), path("${prefix}.collapsed.truncated.fastq.gz"), emit: collapsed_truncated, optional: true
    tuple val(meta), path("${prefix}.paired.fastq.gz")             , emit: paired_interleaved , optional: true
    tuple val(meta), path('*.settings')                            , emit: settings
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def list = adapterlist ? "--adapter-list ${adapterlist}" : ""

    if (meta.single_end) {
        """
        AdapterRemoval  \\
            --file1 ${reads} \\
            ${args} \\
            ${list} \\
            --basename ${prefix} \\
            --threads ${task.cpus} \\
            --seed 42 \\
            --gzip

        ensure_fastq() {
            if [ -f "\${1}" ]; then
                mv "\${1}" "\${1::-3}.fastq.gz"
            fi

        }

        ensure_fastq '${prefix}.truncated.gz'
        ensure_fastq '${prefix}.discarded.gz'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    } else {
        """
        AdapterRemoval  \\
            --file1 ${reads[0]} \\
            --file2 ${reads[1]} \\
            ${args} \\
            ${list} \\
            --basename ${prefix} \\
            --threads ${task.cpus} \\
            --seed 42 \\
            --gzip

        ensure_fastq() {
            if [ -f "\${1}" ]; then
                mv "\${1}" "\${1::-3}.fastq.gz"
            fi

        }

        ensure_fastq '${prefix}.truncated.gz'
        ensure_fastq '${prefix}.discarded.gz'
        ensure_fastq '${prefix}.pair1.truncated.gz'
        ensure_fastq '${prefix}.pair2.truncated.gz'
        ensure_fastq '${prefix}.collapsed.gz'
        ensure_fastq '${prefix}.collapsed.truncated.gz'
        ensure_fastq '${prefix}.paired.gz'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    collapse_cmd = args.contains('--collapse')

    """
    touch '${prefix}.settings'
    echo | gzip > '${prefix}.truncated.fastq.gz'
    echo | gzip > '${prefix}.discarded.fastq.gz'

    if [ "${meta.single_end}" = false ]; then
        echo | gzip > '${prefix}.pair1.truncated.fastq.gz'
        echo | gzip > '${prefix}.pair2.truncated.fastq.gz'
        echo | gzip > '${prefix}.paired.fastq.gz'

        if [ "${collapse_cmd}" = true ]; then
            echo | gzip > '${prefix}.collapsed.truncated.fastq.gz'
            echo | gzip > '${prefix}.collapsed.fastq.gz'
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        adapterremoval: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
    END_VERSIONS
    """
}
