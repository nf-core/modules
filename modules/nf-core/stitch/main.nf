process STITCH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-stitch:1.6.10--r43h06b5641_0':
        'biocontainers/r-stitch:1.6.10--r43h06b5641_0' }"

    input:
    tuple val(meta) , path(posfile), path(input, stageAs: "input"), path(rdata, stageAs: "RData_in"), val(chromosome_name), val(K), val(nGen)
    tuple val(meta2), path(collected_crams), path(collected_crais), path(cramlist)
    tuple val(meta3), path(fasta), path(fasta_fai)
    val seed

    output:
    tuple val(meta), path("input", type: "dir") , emit: input
    tuple val(meta), path("RData", type: "dir") , emit: rdata
    tuple val(meta), path("plots", type: "dir") , emit: plots , optional: { generate_input_only }
    tuple val(meta), path("*.vcf.gz")           , emit: vcf   , optional: { generate_input_only || bgen_output }
    tuple val(meta), path("*.bgen")             , emit: bgen  , optional: { generate_input_only || !bgen_output }
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix               = task.ext.prefix ?: "${meta.id}"
    def args                 = task.ext.args   ?: ""
    def args2                = task.ext.args2  ?: ""
    def generate_input_only  = args2.contains( "--generateInputOnly TRUE" )
    def bgen_output          = args2.contains( "--output_format bgen" )
    def reads_ext            = collected_crams             ? collected_crams.extension.unique()                                : []
    def rsync_cmd            = rdata                       ? "rsync -rL ${rdata}/ RData"                                       : ""
    def stitch_cmd           = seed                        ? "Rscript <(cat \$(which STITCH.R) | tail -n +2 | cat <(echo 'set.seed(${seed})') -)" : "STITCH.R"
    def cramlist_cmd         = cramlist && reads_ext == ["cram"] ? "--cramlist ${cramlist}"                                    : ""
    def bamlist_cmd          = cramlist && reads_ext == ["bam" ] ? "--bamlist ${cramlist}"                                     : ""
    def reference_cmd        = fasta                       ? "--reference ${fasta}"                                            : ""
    def regenerate_input_cmd = input && rdata && !cramlist ? "--regenerateInput FALSE --originalRegionName ${chromosome_name}" : ""
    def rsync_version_cmd    = rdata                       ? "rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //; s/ .*\$//')" : ""
    """
    ${rsync_cmd} ${args}

    ${stitch_cmd} \\
        --chr ${chromosome_name} \\
        --posfile ${posfile} \\
        --outputdir . \\
        --nCores ${task.cpus} \\
        --K ${K} \\
        --nGen ${nGen} \\
        ${cramlist_cmd} \\
        ${bamlist_cmd} \\
        ${reference_cmd} \\
        ${regenerate_input_cmd} \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${rsync_version_cmd}
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-stitch: \$(Rscript -e "cat(as.character(utils::packageVersion('STITCH')))")
    END_VERSIONS
    """

    stub:
    def prefix               = task.ext.prefix      ?: "${meta.id}"
    def args                 = task.ext.args        ?: ""
    def args2                = task.ext.args2       ?: ""
    def generate_input_only  = args2.contains( "--generateInputOnly TRUE" )
    def generate_plots_cmd   = !generate_input_only ? "mkdir plots"                                                                   : ""
    def generate_vcf_cmd     = !generate_input_only ? "touch ${prefix}.vcf.gz"                                                        : ""
    def rsync_version_cmd    = rdata                ? "rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //; s/ .*\$//')" : ""
    """
    touch input
    touch RData
    ${generate_plots_cmd}
    ${generate_vcf_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${rsync_version_cmd}
        r-base: \$(Rscript -e "cat(strsplit(R.version[['version.string']], ' ')[[1]][3])")
        r-stitch: \$(Rscript -e "cat(as.character(utils::packageVersion('STITCH')))")
    END_VERSIONS
    """
}
