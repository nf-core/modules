include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADAPTERREMOVAL {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::adapterremoval=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/adapterremoval:2.3.2--hb7ba0dd_0"
    } else {
        container "quay.io/biocontainers/adapterremoval:2.3.2--hb7ba0dd_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (meta.single_end) {
        """
        AdapterRemoval  \\
            --file1 $reads \\
            $options.args \\
            --basename $prefix \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --output1 ${prefix}.trimmed.fastq.gz \\
            --seed 42 \\
            --gzip \\

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    } else if (!meta.single_end && !meta.collapse) {
        """
        AdapterRemoval  \\
            --file1 ${reads[0]} \\
            --file2 ${reads[1]} \\
            $options.args \\
            --basename $prefix \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --output1 ${prefix}.pair1.trimmed.fastq.gz \\
            --output2 ${prefix}.pair2.trimmed.fastq.gz \\
            --seed 42 \\
            --gzip \\

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    } else {
        """
        AdapterRemoval  \\
            --file1 ${reads[0]} \\
            --file2 ${reads[1]} \\
            --collapse \\
            $options.args \\
            --basename $prefix \\
            --threads $task.cpus \\
            --settings ${prefix}.log \\
            --seed 42 \\
            --gzip \\

        cat *.collapsed.gz *.collapsed.truncated.gz > ${prefix}.merged.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(AdapterRemoval --version 2>&1 | sed -e "s/AdapterRemoval ver. //g")
        END_VERSIONS
        """
    }

}
