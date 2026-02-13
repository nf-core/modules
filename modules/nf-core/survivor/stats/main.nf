process SURVIVOR_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--h9a82719_1':
        'biocontainers/survivor:1.0.7--h9a82719_1' }"

    input:
    tuple val(meta), path(vcf)
    val(minsv)          // Min SV size (-1 to disable)
    val(maxsv)          // Max SV size (-1 to disable)
    val(minnumreads)    // Min number of reads support: RE flag (-1 to disable)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('survivor'), eval("SURVIVOR 2>&1 | grep 'Version' | sed 's/Version: //'"), topic: versions, emit: versions_survivor

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = vcf.getName().endsWith(".gz") ? true : false
    vcf_name = vcf.getName().replace(".gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $vcf > $vcf_name
    fi

    SURVIVOR \\
        stats \\
        $vcf_name \\
        $minsv \\
        $maxsv \\
        $minnumreads \\
        ${prefix}.stats

    rm $vcf_name
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.stats
    """
}
