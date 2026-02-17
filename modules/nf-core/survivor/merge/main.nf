process SURVIVOR_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--h9a82719_1':
        'biocontainers/survivor:1.0.7--h9a82719_1' }"

    input:
    tuple val(meta), path(vcfs)
    val(max_distance_breakpoints)
    val(min_supporting_callers)
    val(account_for_type)
    val(account_for_sv_strands)
    val(estimate_distanced_by_sv_size)
    val(min_sv_size)

    output:
    tuple val(meta), path("*.vcf")   , emit: vcf
    tuple val("${task.process}"), val('survivor'), eval("SURVIVOR 2>&1 | grep 'Version' | sed 's/Version: //'"), topic: versions, emit: versions_survivor

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    vcfs.each{ vcf_file ->
        if ( vcf_file.getExtension() == "gz"){
            error "Gzipped files are not supported by Survivor, please gunzip your VCF files first."
            // https://github.com/fritzsedlazeck/SURVIVOR/issues/158
        }
    }

    """
    SURVIVOR merge \\
        <(ls *.vcf) \\
        ${max_distance_breakpoints} \\
        ${min_supporting_callers} \\
        ${account_for_type} \\
        ${account_for_sv_strands} \\
        ${estimate_distanced_by_sv_size} \\
        ${min_sv_size} \\
        ${prefix}.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    vcfs.each{ vcf_file ->
        if (vcf_file.getExtension() == "gz"){
            error "Gzipped files are not supported by Survivor, please gunzip your VCF files first."
            // https://github.com/fritzsedlazeck/SURVIVOR/issues/158
        }
    }
    """
    touch ${prefix}.vcf
    """
}
