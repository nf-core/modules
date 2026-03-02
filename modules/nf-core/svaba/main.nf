process SVABA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/svaba:1.2.0--h69ac913_1'
        : 'biocontainers/svaba:1.2.0--h69ac913_1'}"

    input:
    tuple val(meta), path(tumorbam), path(tumorbai), path(normalbam), path(normalbai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(bwa_index)
    tuple val(meta5), path(dbsnp)
    tuple val(meta6), path(dbsnp_tbi)
    tuple val(meta7), path(regions)

    output:
    tuple val(meta), path("*.svaba.sv.vcf.gz"), emit: sv, optional: true
    tuple val(meta), path("*.svaba.indel.vcf.gz"), emit: indel, optional: true
    tuple val(meta), path("*.svaba.germline.indel.vcf.gz"), emit: germ_indel, optional: true
    tuple val(meta), path("*.svaba.germline.sv.vcf.gz"), emit: germ_sv, optional: true
    tuple val(meta), path("*.svaba.somatic.indel.vcf.gz"), emit: som_indel, optional: true
    tuple val(meta), path("*.svaba.somatic.sv.vcf.gz"), emit: som_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.sv.vcf.gz"), emit: unfiltered_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.indel.vcf.gz"), emit: unfiltered_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.germline.indel.vcf.gz"), emit: unfiltered_germ_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.germline.sv.vcf.gz"), emit: unfiltered_germ_sv, optional: true
    tuple val(meta), path("*.svaba.unfiltered.somatic.indel.vcf.gz"), emit: unfiltered_som_indel, optional: true
    tuple val(meta), path("*.svaba.unfiltered.somatic.sv.vcf.gz"), emit: unfiltered_som_sv, optional: true
    tuple val(meta), path("*.discordants.txt.gz"), emit: discordants, optional: true
    tuple val(meta), path("*.bps.txt.gz"), emit: raw_calls
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('svaba'), eval("svaba --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n 1"), topic: versions, emit: versions_svaba

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamlist = normalbam ? "-t ${tumorbam} -n ${normalbam}" : "-t ${tumorbam}"
    def dbsnp_cmd = dbsnp ? "--dbsnp-vcf ${dbsnp}" : ""
    def regions_cmd = regions ? "--region ${regions}" : ""
    def bwa = bwa_index ? "cp -s ${bwa_index}/* ." : ""

    """
    ${bwa}

    svaba \\
        run \\
        ${bamlist} \\
        --threads ${task.cpus} \\
        ${dbsnp_cmd} \\
        --id-string ${prefix} \\
        --reference-genome ${fasta} \\
        --g-zip \\
        ${regions_cmd} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.bps.txt.gz
    echo | gzip > ${prefix}.svaba.sv.vcf.gz
    echo | gzip > ${prefix}.svaba.indel.vcf.gz
    echo | gzip > ${prefix}.svaba.germline.indel.vcf.gz
    echo | gzip > ${prefix}.svaba.germline.sv.vcf.gz
    echo | gzip > ${prefix}.svaba.somatic.indel.vcf.gz
    echo | gzip > ${prefix}.svaba.somatic.sv.vcf.gz
    echo | gzip > ${prefix}.unfiltered.sv.vcf.gz
    echo | gzip > ${prefix}.unfiltered.indel.vcf.gz
    echo | gzip > ${prefix}.unfiltered.germline.indel.vcf.gz
    echo | gzip > ${prefix}.unfiltered.germline.sv.vcf.gz
    echo | gzip > ${prefix}.unfiltered.somatic.indel.vcf.gz
    echo | gzip > ${prefix}.unfiltered.somatic.sv.vcf.gz
    echo | gzip > ${prefix}.discondants.txt.gz
    touch ${prefix}.log
    """
}
