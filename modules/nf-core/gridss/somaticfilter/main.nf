process GRIDSS_SOMATICFILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta) , path(vcf)
    tuple val(meta2), path(pondir)

    output:
    tuple val(meta), path("*.high_confidence_somatic.vcf.bgz")    , emit: high_conf_sv
    tuple val(meta), path("*.all_somatic.vcf.bgz")                , emit: all_sv
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pondir_command = pondir ? "--pondir ${pondir}" : ""
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    ## There is a R version issue in the current (2.13.2) bioconda implementation that is hard to pin down
    ## The gridss_somatic_filter wrapper inside the bioconda container hardcodes scriptdir, but we can override that

    #/usr/local/share/gridss-2.13.2-0
    #/opt/conda/share/gridss-2.13.2-4

    #FILE=\$(find / 2>/dev/null -type f -name gridss_somatic_filter)
    #echo \$FILE
    #ORIGDIR=\$(basename \$FILE )
    #echo \$ORIGDIR

    # Find the script and its directory in both conda and docker environments
    SCRIPT_PATH=\$(find / -name libgridss.R 2>/dev/null | head -n 1)

    if [ -z "\$SCRIPT_PATH" ]; then
        echo "Error: libgridss.R not found"
        exit 1
    fi

    # Extract the directory containing the script
    SCRIPT_DIR=\$(dirname "\$SCRIPT_PATH")

    echo "libgridss.R found in directory: \$SCRIPT_DIR"

    # copy the two script files from the container into the working directory
    cp \${SCRIPT_DIR}/libgridss.R .
    cp \${SCRIPT_DIR}/gridss.config.R .

    # end the libgridss.R script based on https://github.com/PapenfussLab/gridss/issues/635
    # This adds as.character() around the expression on line 780.
    sed -i '780s/VariantAnnotation::fixed(vcf)\$ALT/as.character(&)/' libgridss.R

    # need to then pass a different --scriptdir
    \${SCRIPT_DIR}/gridss_somatic_filter \\
        --input $vcf \\
        ${pondir_command} \\
        --output ${prefix}.high_confidence_somatic.vcf \\
        --fulloutput ${prefix}.all_somatic.vcf \\
        --scriptdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "" | bgzip > ${prefix}.high_confidence_somatic.vcf.bgz
    echo "" | bgzip > ${prefix}.all_somatic.vcf.bgz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """
}
