process SURVIVOR_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::survivor=1.0.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--h9a82719_1':
        'quay.io/biocontainers/survivor:1.0.7--h9a82719_1' }"

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.vcf")   , emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: [:]
    def prefix = task.ext.prefix ?: "${meta.id}"

    vcfs.each{
        println(it)
        println(it.getExtension())
        if (it.getExtension() == "gz"){
            error "Gzipped files are not supported by Survivor, please gunzip your VCF files first."
            // https://github.com/fritzsedlazeck/SURVIVOR/issues/158
        }
    }

    """
    ls *.vcf > names.txt

    SURVIVOR merge \\
        names.txt \\
        ${args.max_distance_breakpoints} \\
        ${args.min_supporting_callers} \\
        ${args.account_for_type} \\
        ${args.account_for_sv_strands} \\
        ${args.estimate_distanced_by_sv_size} \\
        ${args.min_sv_size} \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
