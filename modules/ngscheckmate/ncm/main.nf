process NGSCHECKMATE_NCM {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ngscheckmate=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.0--py27r41hdfd78af_0':
        'quay.io/biocontainers/ngscheckmate:1.0.0--py27r41hdfd78af_0' }"

    input:
    path bam
    path snp_bed

    output:
    path "*.pdf"                  , emit: bam
    path "*.txt"                  , emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "output"

    """
    gunzip -f *.vcf.gz

    ls -d "\${PWD}"/*.vcf > vcfList.txt

    ncm.py -V -l vcfList.txt -bed ${snp_bed} -O . -N ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(echo \$(ncm.py --help | grep Identity | sed "s/ *Ensuring Sample Identity v//g"))
    END_VERSIONS
    """
}
