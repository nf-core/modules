process RIBODETECTOR {
	tag "$meta.id"
	label 'process_medium'

	conda "${moduleDir}/environment.yml"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5b/5bd1ef1ec62443f20b84e08612cb008d6abed9647a56422bf71d3174146d8dd1/data' :
        'community.wave.seqera.io/library/ribodetector_pytorch-cuda:8e2cdd88bb757059' }"

	input:
	tuple val(meta), path(fastq)
	val length

	output:
	tuple val(meta), path("*.nonrna*.fastq.gz"), emit: fastq
	tuple val(meta), path("*.log")             , emit: log
	tuple val("${task.process}"), val('ribodetector'), eval('ribodetector --version | sed "s/ribodetector //"'), emit: versions_ribodetector, topic: versions

	when:
	task.ext.when == null || task.ext.when

	script:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"
	ribodetector_bin = task.accelerator ? "ribodetector" : "ribodetector_cpu"
	ribodetector_mem = task.accelerator ? "-m ${task.memory.toGiga()}" : ""
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
	"""

	stub:
	def args = task.ext.args ?: ''
	def prefix = task.ext.prefix ?: "${meta.id}"

	"""
	echo $args

	echo | gzip > ${prefix}.nonrna.1.fastq.gz
	echo | gzip > ${prefix}.nonrna.2.fastq.gz
	touch ${prefix}.log
	"""
}
