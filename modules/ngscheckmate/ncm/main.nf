process NGSCHECKMATE_NCM {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ngscheckmate=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.0--py27r41hdfd78af_0':
        'quay.io/biocontainers/ngscheckmate:1.0.0--py27r41hdfd78af_0' }"

    input:
    path files
    path snp_bed
    val bam_mode

    output:
    path "*.pdf"                  , emit: pdf
    path "*.txt"                  , emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "output"
    def mode_flag = bam_mode ? "-B" : "-V"

    """
    for i in *.vcf.gz; do i2=\$(echo \$i | sed s/.vcf.gz/.vcf/g); gunzip -cdf \$i > \$i2;done

    ncm.py ${mode_flag} -d . -bed ${snp_bed} -O . -N ${prefix} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}
