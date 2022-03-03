process NGSCHECKMATE_NCM {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ngscheckmate=1.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.0--py27r41hdfd78af_1':
        'quay.io/biocontainers/ngscheckmate:1.0.0--py27r41hdfd78af_1' }"

    input:
    path files
    path snp_bed
    path fasta

    output:
    path "*.pdf"                  , emit: pdf
    path "*_corr_matrix.txt"       , emit: corr_matrix
    path "*_matched.txt"           , emit: matched
    path "*_all.txt"           , emit: all
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "output"
    def unzip = false
    opts = args.tokenize()
    if (files.every{ it.toString().endsWith('.bam') || it.toString().endsWith('.bai') } ) {
        if (!opts.contains('-B')) args += ' -B'
        assert !opts.contains('-V')
    } else if ( files.every{ it.toString().endsWith('.vcf') }) {
        if (!opts.contains('-V')) args += ' -V'
        assert !opts.contains('-B')
    } else if ( files.every{ it.toString().endsWith('.vcf.gz') }) {
        unzip = true
        if (!opts.contains('-V')) args += ' -V'
        assert !opts.contains('-B')
    } else {
        throw new Exception("Files must be of the same type")
    }

    """

    if $unzip
    then
        for VCFGZ in *.vcf.gz; do
            gunzip -cdf \$VCFGZ > \$( basename \$VCFGZ .gz );
        done
    fi

    NCM_REF="./"${fasta} ncm.py -d . -bed ${snp_bed} -O . -N ${prefix} $args

    if $unzip
    then
        rm -f *.vcf  # clean up decompressed vcfs
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}
