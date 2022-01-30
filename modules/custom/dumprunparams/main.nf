include { dump_params_yml; indent_code_block } from "./parametrize"

def VERSION="0.0.1"

process CUSTOM_DUMPRUNPARAMS {
    label 'process_low'

    input:
    val(exclude)

    output:
    path "params_mqc.html", emit: mqc_html
    path "versions.yml"  , emit: versions

    script:
    run_params = params

    to_exclude = run_params.subMap(exclude)
    run_params_cleaned = exclude ? params.minus(to_exclude) : params
    run_params_formatted = run_params_cleaned.each{ it.getValue() instanceof LinkedHashMap | it.getValue() instanceof List | it.getValue() instanceof nextflow.config.ConfigMap ? run_params[it.getKey()] = '<truncated_list>' : null }

    yml_cmd = dump_params_yml(run_params_formatted)

    """
    ${yml_cmd}
    echo ""
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dumprunparams: $VERSION
    END_VERSIONS
    """
}
