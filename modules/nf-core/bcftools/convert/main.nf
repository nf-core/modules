process BCFTOOLS_CONVERT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(input), path(input_index)
    tuple val(meta2), path(fasta)
    path(bed)

    output:
    tuple val(meta), path("*.vcf.gz"),      optional:true , emit: vcf_gz
    tuple val(meta), path("*.vcf")   ,      optional:true , emit: vcf
    tuple val(meta), path("*.bcf.gz"),      optional:true , emit: bcf_gz
    tuple val(meta), path("*.bcf")   ,      optional:true , emit: bcf
    tuple val(meta), path("*.hap.gz"),      optional:true , emit: hap
    tuple val(meta), path("*.legend.gz"),   optional:true , emit: legend
    tuple val(meta), path("*.samples"),     optional:true , emit: samples
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def regions = bed ? "--regions-file $bed" : ""
    def reference = fasta ?  "--fasta-ref $fasta" : ""
    def extension = args.contains("--output-type b")    || args.contains("-Ob")    ? "--output ${prefix}.bcf.gz" :
                    args.contains("--output-type u")    || args.contains("-Ou")    ? "--output ${prefix}.bcf" :
                    args.contains("--output-type z")    || args.contains("-Oz")    ? "--output ${prefix}.vcf.gz" :
                    args.contains("--output-type v")    || args.contains("-Ov")    ? "--output ${prefix}.vcf" :
                    args.contains("--haplegendsample")  || args.contains("-h")     ? "" :
                    "--output ${prefix}.vcf.gz"

    """
    bcftools convert \\
        $args \\
        $regions \\
        $extension \\
        --threads $task.cpus \\
        $reference \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf.gz"
    """
    touch ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
