process ONTSPECTRE_CNVCALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/053f18bb8a40ea805d2b1cc594c7b824d3c723d9ba9347dbfcca3832649d0184/data'
        : 'community.wave.seqera.io/library/polars_setuptools_ont-spectre:c268298acfd5e637'}"

    input:
    tuple val(meta), path(summary_txt), path(regions_bed), path(regions_csi), path(snv)
    tuple val(meta2), path(fasta)
    path metadata  // optional
    path blacklist // optional
    val bin_size

    output:
    tuple val(meta), path("${prefix}_cnv.vcf")                , emit: vcf
    tuple val(meta), path("${prefix}_cnv.bed")                , emit: bed
    tuple val(meta), path("${prefix}.spc.gz")                 , emit: spc
    tuple val(meta), path("${prefix}_predicted_karyotype.txt"), emit: karyotype
    tuple val("${task.process}"), val("spectre"), eval("spectre version 2>&1 | tail -n1 | sed 's/.*version: //'"), emit: versions_spectre, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_spectre"
    def metadata_arg = metadata ? "--metadata ${metadata}" : ''
    def blacklist_arg = blacklist ? "--blacklist ${blacklist}" : ''
    """
    mkdir -p coverage
    ln -s ../${summary_txt} coverage/${prefix}.mosdepth.summary.txt
    ln -s ../${regions_bed} coverage/${prefix}.regions.bed.gz
    ln -s ../${regions_csi} coverage/${prefix}.regions.bed.gz.csi

    spectre \\
        CNVCaller \\
        --bin-size ${bin_size} \\
        --coverage coverage/ \\
        --sample-id ${prefix} \\
        --output-dir spectre/ \\
        --reference ${fasta} \\
        --snv ${snv} \\
        --threads $task.cpus \\
        $metadata_arg \\
        $blacklist_arg \\
        $args

    mv spectre/${prefix}.vcf ${prefix}_cnv.vcf
    mv spectre/${prefix}_cnv.bed ${prefix}_cnv.bed
    mv spectre/${prefix}.spc.gz ${prefix}.spc.gz
    mv spectre/predicted_karyotype.txt ${prefix}_predicted_karyotype.txt

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_spectre"
    """
    echo $args
    touch ${prefix}_cnv.vcf
    touch ${prefix}_cnv.bed
    echo "" | gzip > ${prefix}.spc.gz
    touch ${prefix}_predicted_karyotype.txt
    """
}
