process LLAMACPP_PYTHON_RUN {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${task.ext.use_gpu ? 'quay.io/nf-core/llama-cpp-python:0.1.9' : 'community.wave.seqera.io/library/llama-cpp-python:0.3.16--b351398cd0ea7fc5'}"

    input:
    tuple val(meta), path(prompt_file), path(gguf_model)

    output:
    tuple val(meta), path("${prefix}.txt"), emit: output
    tuple val("${task.process}"), val("llama-cpp-python"), eval("python3 -c 'import llama_cpp; print(llama_cpp.__version__)'"), topic: versions, emit: versions_llama_cpp_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    llamacpp-python.py \
        --model ${gguf_model} \
        --messages ${prompt_file} \
        --output ${prefix}.txt \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
