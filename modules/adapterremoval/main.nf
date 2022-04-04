process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::adapterremoval=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/adapterremoval:2.3.2--hb7ba0dd_0' :
        'quay.io/biocontainers/adapterremoval:2.3.2--hb7ba0dd_0' }"

    input:
    tuple val(meta), path(reads)
    path(adapterlist)

    output:
    tuple val(meta), path("${prefix}.truncated.fastq.gz")            , optional: true, emit: singles_truncated
    tuple val(meta), path("${prefix}.discarded.fastq.gz")            , optional: true, emit: discarded
    tuple val(meta), path("${prefix}.pair{1,2}.truncated.fastq.gz")  , optional: true, emit: paired_truncated
    tuple val(meta), path("${prefix}.collapsed.fastq.gz")            , optional: true, emit: collapsed
    tuple val(meta), path("${prefix}.collapsed.truncated.fastq.gz")  , optional: true, emit: collapsed_truncated
    tuple val(meta), path("${prefix}.paired.fastq.gz")               , optional: true, emit: paired_interleaved
    tuple val(meta), path('*.settings')                              , emit: settings
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def list = adapterlist ? "--adapter-list ${adapterlist}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        AdapterRemoval  \\
            --file1 $reads \\
            $args \\
            $adapterlist \\
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
            $args \\
            $adapterlist \\
            --basename ${prefix} \\
            --threads $task.cpus \\
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

}
