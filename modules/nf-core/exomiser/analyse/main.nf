process EXOMISER_ANALYSE {
    tag "${meta.id}"
    label 'process_medium'

    container "nf-core/exomiser-cli:15.0.0-bash"

    input:
    tuple val(meta), path(vcf), path(ped), val(assembly), path(phenopacket), path(analysis_script)
    tuple val(meta2), path(reference_cache, stageAs: 'exomiser_data/*'), val(reference_version)
    tuple val(meta3), path(phenotype_cache, stageAs: 'exomiser_data/*'), val(phenotype_version)

    output:
    tuple val(meta), path("*.{tsv}"), emit: tsv
    tuple val(meta), path("*.{json}"), emit: json
    tuple val(meta), path("*.{html}"), emit: html
    tuple val(meta), path("*.{parquet}"), emit: parquet
    tuple val(meta), path("*.{vcf}"), emit: vcf
    tuple val("${task.process}"), val('exomiser'), eval("exomiser --version"), topic: versions, emit: versions_exomiser

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("EXOMISERCLI_ANALYSE module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ped_cmd = ped ? "--ped=${ped}" : ""
    def phenopacket_cmd = phenopacket ? "--sample=${phenopacket}" : ""
    def assembly_cmd = assembly ? "--assembly=${assembly}" : ""
    def analysis_cmd = analysis_script ? "--analysis ${analysis_script}" : ""
    def vcf_cmd = vcf ? "--vcf=${vcf}" : ""

    """
        export EXOMISER_DATA_DIRECTORY=./exomiser_data
        export EXOMISER_${assembly}_DATA_VERSION=${reference_version}
        export EXOMISER_PHENOTYPE_DATA_VERSION=${phenotype_version}

        exomiser analyse \\
        ${ped_cmd} \\
        ${phenopacket_cmd} \\
        ${assembly_cmd} \\
        ${vcf_cmd}\\
        ${analysis_cmd} \\
        ${args} \\
        --output-directory=\$PWD \
        --output-filename=${prefix} \
        --exomiser.data-directory=./exomiser-data \
        --exomiser.${assembly}.data-version=${reference_version} \
        --exomiser.phenotype.data-version=${phenotype_version}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}
    touch ${prefix}.tsv
    touch ${prefix}.json
    touch ${prefix}.html
    touch ${prefix}.parquet
    touch ${prefix}.vcf
    """
}
