process LLAMACPPPYTHON_RUN {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${task.accelerator ? 'quay.io/nf-core/llama-cpp-python:0.1.9' : (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'oras://community.wave.seqera.io/library/llama-cpp-python:0.3.16--d6f959a4c13960c4' : 'community.wave.seqera.io/library/llama-cpp-python:0.3.16--b351398cd0ea7fc5')}"

    input:
    tuple val(meta), path(prompt_file), path(gguf_model)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: output
    path "versions.yml", emit: versions_llama_cpp_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    template('llama-cpp-python.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        llama-cpp-python: \$(python3 -c 'import llama_cpp; print(llama_cpp.__version__)')
    END_VERSIONS
    """
}
