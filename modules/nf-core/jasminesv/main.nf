process JASMINESV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs), path(bams), path(sample_dists)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(chr_norm)

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def args3   = task.ext.args3 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    def make_bam = bams ? "ls *.bam > bams.txt" : ""
    def bam_argument = bams ? "bam_list=bams.txt" : ""
    def iris_argument = args2 != '' ? "iris_args=${args2}" : ""
    def sample_dists_argument = sample_dists ? "sample_dists=${sample_dists}" : ""
    def chr_norm_argument = chr_norm ? "chr_norm_file=${chr_norm}" : ""

    vcfs.each { vcf ->
        if ("$vcf".startsWith("${prefix}.vcf")) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    def unzip_inputs = vcfs.collect { it.baseName.endsWith(".vcf.gz") ? "    bgzip -d --threads ${task.cpus} ${args2} ${it}" : "" }.join("\n")
    """
    ${unzip_inputs}

    ls *.vcf > vcfs.txt
    ${make_bam}

    jasmine \\
        file_list=vcfs.txt \\
        out_file=${prefix}.vcf \\
        threads=${task.cpus} \\
        genome_file=${fasta} \\
        ${bam_argument} \\
        ${iris_argument} \\
        ${sample_dists_argument} \\
        ${chr_norm_argument} \\
        ${args}

    bgzip --threads ${task.cpus} ${args3} ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(jasmine 2>&1 | grep "version" | sed 's/Jasmine version //'))
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(jasmine 2>&1 | grep "version" | sed 's/Jasmine version //'))
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
