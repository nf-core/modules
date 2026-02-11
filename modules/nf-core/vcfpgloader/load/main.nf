// Database connection parameters are passed as val() inputs to support
// dynamic per-sample configuration. Set PGPASSWORD via environment variable
// or Nextflow secrets before running:
//
// Option 1 - Environment variable:
//   export PGPASSWORD='your_password'
//
// Option 2 - Nextflow secrets (nextflow.config):
//   env {
//       PGPASSWORD = secrets.PGPASSWORD
//   }

process VCFPGLOADER_LOAD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf-pg-loader:0.5.4--pyhdfd78af_0' :
        'biocontainers/vcf-pg-loader:0.5.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi), val(db_host), val(db_port), val(db_name), val(db_user), val(db_schema)

    output:
    tuple val(meta), path("*.load_report.json"), emit: report
    tuple val(meta), path("*.load.log"), emit: log
    tuple val(meta), env('ROWS_LOADED'), emit: row_count
    tuple val("${task.process}"), val("vcf-pg-loader"), eval("vcf-pg-loader --version | sed 's/.*version //'"), topic: versions, emit: versions_vcfpgloader

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // NOTE: batch_size exposed via task.ext for pipeline-level tuning of memory/performance tradeoffs
    def batch_size = task.ext.batch_size ?: '10000'
    """
    vcf-pg-loader load \\
        --host ${db_host} \\
        --port ${db_port} \\
        --database ${db_name} \\
        --user ${db_user} \\
        --schema ${db_schema} \\
        --batch ${batch_size} \\
        --workers ${task.cpus} \\
        --sample-id ${meta.id} \\
        --report ${prefix}.load_report.json \\
        --log ${prefix}.load.log \\
        ${args} \\
        ${vcf}

    export ROWS_LOADED=\$(python3 -c "import json; print(json.load(open('${prefix}.load_report.json'))['variants_loaded'])")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_JSON > ${prefix}.load_report.json
    {"status": "stub", "variants_loaded": 0, "elapsed_seconds": 0}
    END_JSON
    touch ${prefix}.load.log
    export ROWS_LOADED=0
    """
}
