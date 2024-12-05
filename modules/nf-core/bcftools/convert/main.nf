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
    tuple val(meta), path("*.samples")  ,   optional:true , emit: samples
    tuple val(meta), path("*.tbi")      ,   optional:true , emit: tbi
    tuple val(meta), path("*.csi")      ,   optional:true , emit: csi
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def regions   = bed ? "--regions-file $bed" : ""
    def reference = fasta ?  "--fasta-ref $fasta" : ""
    def extension = args.contains("--output-type b")   || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u")   || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z")   || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v")   || args.contains("-Ov") ? "vcf" :
                    args.contains("--haplegendsample") || args.contains("-h")  ? "" :
                    "vcf.gz"

    def output_cmd = args.contains("--haplegendsample") ? "" : "--output ${prefix}.${extension}"

    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    bcftools convert \\
        $args \\
        $regions \\
        $output_cmd \\
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

    def extension = args.contains("--output-type b")   || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u")   || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z")   || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v")   || args.contains("-Ov") ? "vcf" :
                    args.contains("--haplegendsample") || args.contains("-h")  ? "hap.gz" :
                    "vcf.gz"

    def index = args.contains("--write-index=tbi") || args.contains("-W=tbi") ? "tbi" :
            args.contains("--write-index=csi") || args.contains("-W=csi") ? "csi" :
            args.contains("--write-index") || args.contains("-W") ? "csi" :
            ""

    def create_cmd   = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    def hap = ""

    if (args.contains('--haplegendsample')) {
        def args_split = args.split(' ')
        hap = args_split.findIndexOf{ it == '--haplegendsample'}
        prefix = args_split[hap + 1]
    }
    """
    if [ -n "${hap}" ] ; then
        ${create_cmd} ${prefix}.hap.gz
        ${create_cmd} ${prefix}.legend.gz
        touch ${prefix}.samples
    else
        ${create_cmd} ${prefix}.${extension}
        ${create_index}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
