process JASMINESV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs, arity:'1..*'), path(bams), path(bais), path(sample_dists)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    path(chr_norm)

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val("${task.process}"), val("jasminesv"), eval('jasmine 2>&1 | grep "version" | sed "s/Jasmine version //"'), topic: versions, emit: versions_jasminesv
    tuple val("${task.process}"), val("bgzip"), eval('bgzip --version | head -1 | sed -e "s/bgzip (htslib) //"'), topic: versions, emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def args3   = task.ext.args3 ?: ''
    def args4   = task.ext.args4 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    def make_bam = bams ? "ls *.bam > bams.txt" : ""
    def bam_argument = bams ? "bam_list=bams.txt" : ""
    def iris_argument = args2 != '' ? "iris_args=${args2}" : ""
    def sample_dists_argument = sample_dists ? "sample_dists=${sample_dists}" : ""
    def chr_norm_argument = chr_norm ? "chr_norm_file=${chr_norm}" : ""

    def first_vcf = vcfs[0].name.replaceAll(".gz\$", "")
    def unzip_inputs = vcfs.collect { vcf -> vcf.name.endsWith(".vcf.gz") ? "    bgzip -d --threads ${task.cpus} ${args4} ${vcf}" : "" }.join("\n")

    vcfs.each { vcf ->
        if ("$vcf".startsWith("${prefix}.vcf")) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

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

    if [ -s ${prefix}.vcf ]; then
        echo "The file is not empty"
        bgzip --threads ${task.cpus} ${args3} ${prefix}.vcf
    else
        echo "The file is empty, using the header of the first VCF as output file"
        cat ${first_vcf} | grep "#" | bgzip --threads ${task.cpus} ${args3} > ${prefix}.vcf.gz
    fi
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    vcfs.each { vcf ->
        if ("$vcf".startsWith("${prefix}.vcf")) error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
