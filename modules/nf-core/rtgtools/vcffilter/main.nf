process RTGTOOLS_VCFFILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/rtg-tools:3.12.1--hdfd78af_0'
        : 'biocontainers/rtg-tools:3.12.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path bed
    path vcffilter_js
    path optional_files

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), emit: filtered_vcf
    tuple val(meta), path("${meta.id}.filtered.vcf.gz.tbi"), emit: filtered_index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def js_filter = vcffilter_js

    def filter_func_STR = task.ext.filter_func_STR ? "--keep-expr '${task.ext.filter_func_STR}'" : ""
    def filter_func_JS_path = vcffilter_js ? "--javascript ${vcffilter_js}" : ""
    def filter_func_arg = filter_func_STR ?: filter_func_JS_path

    def bed_regions = bed ? "--bed-regions ${bed}" : ""
    def index_vcf = (bed && !tbi) ? "rtg index ${vcf}" : ""

    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_vcf = "${prefix}_${task.ext.output_name}_filtered.vcf.gz"

    """
    rtg vcffilter -i ${vcf} -o ${output_vcf} ${args} ${filter_func_arg} ${bed_regions}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rtg-tools: \$(echo \$(rtg version | head -n 1 | awk '{print \$4}'))
    END_VERSIONS
    """
}

