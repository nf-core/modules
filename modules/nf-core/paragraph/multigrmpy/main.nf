process PARAGRAPH_MULTIGRMPY {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::paragraph=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paragraph:2.3--h21f15d8_1':
        'biocontainers/paragraph:2.3--h21f15d8_1' }"

    input:
    tuple val(meta), path(variants), path(variants_index), path(reads), path(reads_index), path(manifest)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    tuple val(meta), path("*.json.gz")  , emit: json, optional:true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    def check_vcf = variants.name.endsWith(".vcf.gz") ? "variant=\$(bgzip -d --threads ${task.cpus} --stdout ${variants} | awk '/^#/ {next} {print 1;exit}' || echo 0)":
                    variants.extension == "vcf" ? "variant=\$(cat ${variants} | awk '/^#/ {next} {print 1;exit}' || echo 0)":
                    "variant=1"

    if ("${variants}" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("${variants}" == "${prefix}.json.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    ${check_vcf}

    if [ \$variant -eq 1 ]
    then
        multigrmpy.py \\
            --input ${variants} \\
            --manifest ${manifest} \\
            --output ${prefix} \\
            --reference ${fasta} \\
            --threads ${task.cpus} \\
            ${args}

        mv ${prefix}/genotypes.vcf.gz ${prefix}.vcf.gz
        mv ${prefix}/genotypes.json.gz ${prefix}.json.gz
    else
        echo "${variants} was empty, so the multigrmpy.py process was skipped."
        cp ${variants} ${prefix}.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    if ("${variants}" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    if ("${variants}" == "${prefix}.json.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.json.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paragraph: ${VERSION}
    END_VERSIONS
    """
}
