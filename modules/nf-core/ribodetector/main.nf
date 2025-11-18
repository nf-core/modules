process RIBODETECTOR {
	tag "$meta.id"
	label 'process_medium'

	conda "${moduleDir}/environment.yml"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
		'https://depot.galaxyproject.org/singularity/ribodetector:0.3.1--pyhdfd78af_0':
		'biocontainers/ribodetector:0.3.1--pyhdfd78af_0' }"

	input:
	tuple val(meta), path(fastq)
	val length

	output:
	tuple val(meta), path("*.nonrna*.fastq.gz"), emit: fastq
	tuple val(meta), path("*.log")             , emit: log
	path "versions.yml"                        , emit: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"
	ribodetector_bin = task.accelerator ? "ribodetector" : "ribodetector_cpu"
	ribodetector_mem = task.accelerator ? "-m $task.memory.toGiga()" : ""
	output = meta.single_end ? "${prefix}.nonrna.fastq.gz" : "${prefix}.nonrna.1.fastq.gz ${prefix}.nonrna.2.fastq.gz"

	"""
	${ribodetector_bin} \\
		-i ${fastq} \\
		-o ${output} \\
		-l ${length} \\
		-t ${task.cpus} \\
		--log ${prefix}.log \\
		${ribodetector_mem} \\
		${args}

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		ribodetector: \$(ribodetector --version | sed 's/ribodetector //g')
	END_VERSIONS
	"""

	stub:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"

	"""
	echo $args

	echo | gzip > ${prefix}.nonrna.1.fastq.gz
    echo | gzip > ${prefix}.nonrna.2.fastq.gz
	touch ${prefix}.log

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
		ribodetector: \$(ribodetector --version | sed 's/ribodetector //g')
	END_VERSIONS
	"""
}
