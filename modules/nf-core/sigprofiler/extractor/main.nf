process SIGPROFILER_EXTRACTOR {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: List required Conda package(s). Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(vcf)
    val(genome)

    output:
    tuple val(meta), path("results"), emit: results
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "results"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    """
    #!/usr/bin/env python
    from SigProfilerExtractor import sigpro as sig
    from SigProfilerMatrixGenerator import install as genInstall

    genInstall.install(${genome})
    sig.sigProfilerExtractor("vcf", ${outdir}, ${vcf}, minimum_signatures=1, maximum_signatures=3)

    # Versions
    f = open("versions.yml", "a")
    f.write('"${task.process}":\\n')
    f.write("    SigProfilerExtractor: " + SigProfilerExtractor.__version__ + "\\n")
    f.write("    SigProfilerMatrixGenerator: " + SigProfilerMatrixGenerator.__version__ + "\\n")
    f.close()
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outdir = "results"
    // Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    // Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sigprofiler: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
