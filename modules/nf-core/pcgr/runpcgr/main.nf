process PCGR {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pcgr:2.2.1--6e74968e17eb9f97':
        'community.wave.seqera.io/library/pcgr:2.2.1--cf53926ac45a4bda' }"

    input:
    tuple val(meta), path(vcf) // Maybe tbi is needed for vcf
    path (vep_dir)
    path (refdata_dir)
    path (cna)
    path (cpsr)
    path (cpsr_yaml)

    output:
    tuple val(meta), path ("${meta.id}"), emit: pcgr_reports
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def sample_id = task.ext.sample_id ?: "$prefix"
    def cna_option    = cna ? "--input_cna $cna" : ''
    def cpsr_option   = cpsr ? "--input_cpsr $cpsr" : ''
    def cpsr_yaml_option = cpsr_yaml ? "--input_cpsr_yaml $cpsr_yaml" : ''

    """
    mkdir -p $prefix

    pcgr \\
        --input_vcf $vcf \\
        --vep_dir $vep_dir \\
        --refdata_dir $refdata_dir \\
        --output_dir $prefix \\
        --genome_assembly $genome \\
        --sample_id $sample_id \\
        $cna_option \\
        $cpsr_option \\
        $cpsr_yaml_option \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pcgr: \$(echo \$( pcgr --version | sed 's/pcgr //g' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // if (cna) {
    //     """
    //     touch ${prefix}/${prefix}.cna_gene.tsv.gz
    //     touch ${prefix}/${prefix}.cna_gene_ann.tsv.gz
    //     touch ${prefix}/${prefix}.cna_segment.tsv.gz
    //     """
    // }
    // if (args?.contains('--vcf2maf')) {
    //     """
    //     touch ${prefix}/${prefix}.maf
    //     """
    // }
    // if (!args?.contains('--no_html')) {
    //     """
    //     touch ${prefix}/${prefix}.html
    //     """
    // }

    // """

    """
    echo $args
    mkdir -p $prefix

    touch ${prefix}/${prefix}.conf.yaml
    touch ${prefix}/${prefix}.expression.tsv.gz
    touch ${prefix}/${prefix}.expression_outliers.tsv.gz
    touch ${prefix}/${prefix}.expression_similarity.tsv.gz
    touch ${prefix}/${prefix}.log
    touch ${prefix}/${prefix}.msigs.tsv.gz
    touch ${prefix}/${prefix}.pass.tsv.gz
    touch ${prefix}/${prefix}.pass.vcf.gz
    touch ${prefix}/${prefix}.snv_indel_ann.tsv.gz
    touch ${prefix}/${prefix}.pass.vcf.gz.tbi
    touch ${prefix}/${prefix}.tmb.tsv
    touch ${prefix}/${prefix}.vcf.gz
    touch ${prefix}/${prefix}.vcf.gz.tbi
    touch ${prefix}/${prefix}.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    pcgrtemplate: \$(pcgrtemplate --version)
    END_VERSIONS
    """
}
