process PCGR {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pcgr:2.2.1--6e74968e17eb9f97':
        'community.wave.seqera.io/library/pcgr:2.2.1--cf53926ac45a4bda' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna)
    path(pcgr_dir), stageAs: "PCGR/data/${params.genome.toLowerCase()}"
    path vep_cache

    output:
    tuple val(meta), path("${prefix}"), emit: pcgr_reports
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def database = './PCGR'
    def args     = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def cna_param      = params.cna_analysis ? "--input_cna $cna" : ''
    """
    export XDG_CACHE_HOME=/tmp
    export XDG_DATA_HOME=/tmp
    export QUARTO_PRINT_STACK=true

    mkdir -p $prefix

    pcgr \\
        --input_vcf $vcf \\
        --vep_dir $vep_cache \\
        --refdata_dir $database \\
        --output_dir $prefix \\
        --genome_assembly $genome \\
        --sample_id $prefix \\
        --tumor_dp_tag 'TDP' \\
        --tumor_af_tag 'TAF' \\
        --call_conf_tag 'TAL' \\
        $cna_param \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr //g' ))
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args
    
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgrtemplate: \$(pcgrtemplate --version)
    END_VERSIONS
    """
}
