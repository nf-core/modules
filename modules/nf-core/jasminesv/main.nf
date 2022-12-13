process JASMINESV {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::jasminesv=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs), path(bams), path(sample_dists)
    path(fasta)
    path(chr_norm)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def args2   = task.ext.args2 ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    vcfs.each{
        println(it)
        println(it.getExtension())
        if (it.getExtension() == "gz"){
            error "Gzipped files are not supported by Jasmine, please gunzip your VCF files first."
            // https://github.com/mkirsche/Jasmine/issues/31
        }
    }

    make_bam = bams ? "ls *.bam > bams.txt" : ""
    bam_argument = bams ? "bam_list=bams.txt"
    make_iris = args2 != '' ? "echo ${args2} > iris.txt"
    iris_argument = args2 != '' ? "iris_args=iris.txt"
    sample_dists_argument = sample_dists ? "sample_dists=${sample_dists}" : ""
    chr_norm_argument = chr_norm ? "chr_norm_file=${chr_norm}" : ""

    """
    ls *.vcf > vcfs.txt
    ${make_bam}
    ${make_iris}

    jasmine \\
        file_list=vcfs.txt \\
        out_file=${prefix}.vcf \\
        threads=${task.cpus} \\
        genome_file=${fasta} \\
        ${bam_argument} \\
        ${iris_argument} \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
