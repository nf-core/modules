// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FREEBAYES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3"
    } else {
        container "quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3"
    }

    input:
    tuple val(meta), path(bam)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (task.cpus > 0) {
        """
        freebayes-parallel \\
            <(fasta_generate_regions.py ${fasta}.fai 10000) ${task.cpus} \
            -f $fasta \\
            $options.args \\
            $bam  > ${prefix}.vcf

        gzip ${prefix}.vcf
        echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' > ${software}.version.txt
        """

    } else {
        """
        freebayes \\
            -f $fasta \\
            $options.args \\
            $bam > ${prefix}.vcf

        gzip ${prefix}.vcf
        echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' > ${software}.version.txt
        """
    }
}
