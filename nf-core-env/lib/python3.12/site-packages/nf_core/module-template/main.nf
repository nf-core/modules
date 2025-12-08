{%- if not_empty_template -%}
// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.
{%- endif %}

process {{ component_name_underscore|upper }} {
    tag {{ '"$meta.id"' if has_meta else "'$bam'" }}
    label '{{ process_label }}'

    {% if not_empty_template -%}
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    {% endif -%}
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '{{ singularity_container if singularity_container else 'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE' }}':
        '{{ docker_container if docker_container else 'biocontainers/YOUR-TOOL-HERE' }}' }"

    input:
    {%- if inputs %}
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    {%- for input_name, ontologies in inputs.items() %}
    {% set meta_index = loop.index|string if not loop.first else '' %}
    {{ 'tuple val(meta' + meta_index + '), path(' + input_name + ')' if has_meta else 'path ' + input_name }}
    {%- endfor %}
    {%- else -%}
    {% if not_empty_template -%}
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    {{ 'tuple val(meta), path(bam)' if has_meta else 'path bam' }}
    {%- else -%}
    {{ 'tuple val(meta), path(input)' if has_meta else 'path input' }}
    {%- endif %}
    {%- endif %}

    output:
    {%- if outputs %}
    // TODO nf-core: Update the information obtained from bio.tools and make sure that it is correct
    {%- for output_name, ontologies in outputs.items() %}
    {{ 'tuple val(meta), path("*.{' + ontologies[2]|join(',') + '}")' if has_meta else 'path ' + output_name }}, emit: {{ output_name }}
    {%- endfor %}
    {%- else %}
    {% if not_empty_template -%}
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    {{ 'tuple val(meta), path("*.bam")' if has_meta else 'path "*.bam"' }}, emit: bam
    // TODO nf-core: List additional required output channels/values here
    {%- else -%}
    {{ 'tuple val(meta), path("*")' if has_meta else 'path "*"' }}, emit: output
    {%- endif %}
    {%- endif %}
    {% if not_empty_template -%}
    // TODO nf-core: Update the command here to obtain the version number of the software used in this module
    // TODO nf-core: If multiple software packages are used in this module, all MUST be added here
    //               by copying the line below and replacing the current tool with the extra tool(s)
    {%- endif %}
    tuple val("${task.process}"), val('{{ component }}'), eval("{{ component }} --version"), topic: versions, emit: versions_{{ component }}

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    {% if has_meta -%}
    def prefix = task.ext.prefix ?: "${meta.id}"
    {%- endif %}
    {% if not_empty_template -%}
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    {%- endif %}
    """
    {% if not_empty_template -%}
    {{ component }} \\
        $args \\
        -@ $task.cpus \\
        {%- if has_meta %}
        {%- if inputs %}
        {%- for input_name, ontologies in inputs.items() %}
        {%- set extensions = ontologies[2] %}
        {%- for ext in extensions %}
        -o ${prefix}.{{ ext }} \\
        {%- endfor %}
        {%- endfor %}
        {%- else %}
        -o ${prefix}.bam \\
        {%- endif %}
        {%- endif %}
        {%- if inputs %}
        {%- for input_name, ontologies in inputs.items() %}
        ${{ input_name }} \\
        {%- endfor %}
        {%- else %}
        $bam
        {%- endif %}
    {%- endif %}
    """

    stub:
    def args = task.ext.args ?: ''
    {% if has_meta -%}
    def prefix = task.ext.prefix ?: "${meta.id}"
    {%- endif %}
    {% if not_empty_template -%}
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    {%- endif %}
    """
    echo $args
    {% if not_empty_template -%}
    {%- if inputs %}
    {%- for input_name, ontologies in inputs.items() %}
    {%- set extensions = ontologies[2] %}
    {%- for ext in extensions %}
    touch ${prefix}.{{ ext }}
    {%- endfor %}
    {%- endfor %}
    {%- else %}
    touch ${prefix}.bam
    {%- endif %}
    {%- endif %}
    """
}
