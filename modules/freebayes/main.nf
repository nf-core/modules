// Import generic module functions
include { initOptions; saveFiles; getProcessName; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FREEBAYES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3"
    } else {
        container "quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3"
    }

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (task.cpus > 0) {
        """
        freebayes-parallel \\
            <(fasta_generate_regions.py ${fasta}.fai 10000) ${task.cpus} \
            -f $fasta \\
            $options.args \\
            $bam  > ${prefix}.vcf

        gzip ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
        END_VERSIONS
        """

    } else {
        """
        freebayes \\
            -f $fasta \\
            $options.args \\
            $bam > ${prefix}.vcf

        gzip ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
        END_VERSIONS
        """
    }
}
