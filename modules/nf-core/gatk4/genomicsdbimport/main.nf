process GATK4_GENOMICSDBIMPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e1/e1330cb37b3acde07db0a04befa76973fea72d1d10a79542132e3d77950225af/data'
        : 'community.wave.seqera.io/library/gatk4-main_gcnvkernel:b945148cc5e2a616'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(interval_file), val(interval_value), path(wspace)
    val run_intlist
    val run_updatewspace
    val input_map

    output:
    tuple val(meta), path("${prefix}"), emit: genomicsdb, optional: true
    tuple val(meta), path("${updated_db}"), emit: updatedb, optional: true
    tuple val(meta), path("*.interval_list"), emit: intervallist, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // settings for running default create gendb mode
    input_command = input_map ? "--sample-name-map ${vcf[0]}" : vcf.collect { vcf_ -> "--variant ${vcf_}" }.join(' ')

    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    interval_command = interval_file ? "--intervals ${interval_file}" : "--intervals ${interval_value}"
    updated_db = ""

    // settings changed for running get intervals list mode if run_intlist is true
    if (run_intlist) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    }

    // settings changed for running update gendb mode. input_command same as default, update_db forces module to emit the updated gendb
    if (run_updatewspace) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = ''
        updated_db = "${wspace}"
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK GenomicsDBImport] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        GenomicsDBImport \\
        ${input_command} \\
        ${genomicsdb_command} \\
        ${interval_command} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    interval_command = interval_file ? "--intervals ${interval_file}" : "--intervals ${interval_value}"
    updated_db = ""

    // settings changed for running get intervals list mode if run_intlist is true
    if (run_intlist) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    }

    // settings changed for running update gendb mode. input_command same as default, update_db forces module to emit the updated gendb
    if (run_updatewspace) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = ''
        updated_db = "${wspace}"
    }

    def stub_genomicsdb = genomicsdb_command == "--genomicsdb-workspace-path ${prefix}" ? "touch ${prefix}" : ""
    def stub_interval = interval_command == "--output-interval-list-to-file ${prefix}.interval_list" ? "touch ${prefix}.interval_list" : ""
    def stub_update = updated_db != "" ? "touch ${wspace}" : ""

    """
    ${stub_genomicsdb}
    ${stub_interval}
    ${stub_update}
    """
}
