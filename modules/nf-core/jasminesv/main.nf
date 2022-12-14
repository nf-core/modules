process JASMINESV {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::jasminesv=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jasminesv:1.1.5--hdfd78af_0':
        'quay.io/biocontainers/jasminesv:1.1.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcfs), path(bams), path(sample_dists)
    path(fasta)
    path(fasta_fai)
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
        if (it.getExtension() == "gz"){
            error "Gzipped files are not supported by Jasmine, please gunzip your VCF files first."
            // https://github.com/mkirsche/Jasmine/issues/31
        }
    }

    make_bam = bams ? "ls *.bam > bams.txt" : ""
    bam_argument = bams ? "bam_list=bams.txt" : ""
    iris_argument = args2 != '' ? "iris_args=${args2}" : ""
    sample_dists_argument = sample_dists ? "sample_dists=${sample_dists}" : ""
    chr_norm_argument = chr_norm ? "chr_norm_file=${chr_norm}" : ""

    """
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


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(jasmine 2>&1 | grep "version" | sed 's/Jasmine version //'))
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"

    vcfs.each{
        if (it.getExtension() == "gz"){
            error "Gzipped files are not supported by Jasmine, please gunzip your VCF files first."
            // https://github.com/mkirsche/Jasmine/issues/31
        }
    }

    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jasminesv: \$(echo \$(jasmine 2>&1 | grep "version" | sed 's/Jasmine version //'))
    END_VERSIONS
    """
}
