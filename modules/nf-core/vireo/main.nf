// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VIREO {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2eac294d6de9fdecf7b2bede5798788af111ac10ad6f3abc72a45fb7f317af92/data'
:         'community.wave.seqera.io/library/vireosnp:0.5.8--376a9219dbddd077' }"

    input:
    tuple val(meta), path(cell_data), val(n_donor), path(donor_file), path(vartrix_data)
    output:
    tuple val(meta), path('*_summary.tsv')           , emit: summary
    tuple val(meta), path('*_donor_ids.tsv')         , emit: donor_ids
    tuple val(meta), path('*_prob_singlet.tsv.gz')   , emit: prob_singlets
    tuple val(meta), path('*_prob_doublet.tsv.gz')   , emit: prob_doublets
    tuple val(meta), path('*_GT_donors.vireo.vcf.gz'), emit: genotype_vcf     , optional: true
    tuple val(meta), path('*_filtered_variants.tsv') , emit: filtered_variants, optional: true
    tuple val("${task.process}"), val('vireo'), eval('vireo | sed "1!d ; s/Welcome to vireoSNP v//; s/\\!//"'), emit: versions_vireo, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    // hardcode a default random seed
    if (!(args ==~ /.*--randSeed.*/)) {args += " --randSeed 42"}
    // use the same randSeed of vireo for GTbarcode if specified in args
    def matcher = (args =~ /(--randSeed\s+\d+)/)
    def randSeed_GTbarcode = matcher ? matcher[0][1] : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = cell_data       ? "-c ${cell_data}" : "--vartrixData ${vartrix_data}"

    """
    vireo \\
        ${input} \\
        -N ${n_donor} \\
        -d ${donor_file} \\
        -p $task.cpus \\
        -o . \\
        $args

    mv summary.tsv ${prefix}_summary.tsv
    mv donor_ids.tsv ${prefix}_donor_ids.tsv
    mv prob_singlet.tsv.gz ${prefix}_prob_singlet.tsv.gz
    mv prob_doublet.tsv.gz ${prefix}_prob_doublet.tsv.gz
    if [[ -f GT_donors.vireo.vcf.gz ]]; then
        mv GT_donors.vireo.vcf.gz "${prefix}_GT_donors.vireo.vcf.gz"
        GTbarcode -i ${prefix}_GT_donors.vireo.vcf.gz -o ./${prefix}_filtered_variants.tsv ${randSeed_GTbarcode}
    fi
    """

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def optional_files = ''
    if (args.contains('--forceLearnGT')) {
        optional_files = """
        echo "" | gzip > ${prefix}_GT_donors.vireo.vcf.gz
        touch ${prefix}_filtered_variants.tsv
        """
    }

    """
    touch ${prefix}_summary.tsv
    touch ${prefix}_donor_ids.tsv
    echo "" | gzip > ${prefix}_prob_singlet.tsv.gz
    echo "" | gzip > ${prefix}_prob_doublet.tsv.gz

    ${optional_files}
    """
}
