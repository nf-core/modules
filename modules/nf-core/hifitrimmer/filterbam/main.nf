// Adapte from https://github.com/nf-core/modules/tree/master/modules/nf-core/pb/hifitrimmer/filterbam

process HIFITRIMMER_FILTERBAM {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/ddcbb21fe5a0ddaaf0fdccedc5c60d6298710e87eea39a45a75969528f001e37/data' :
      'community.wave.seqera.io/library/hifi_trimmer_htslib_samtools:faf8da27d58784cd' }"

   input:
   tuple val(meta), path(input), path(bed)


   output:
   tuple val(meta), path("*.fast{q,a}.gz"), emit: filtered
   tuple val("${task.process}"), val('hifi_trimmer'), eval("hifi_trimmer --version | cut -d' ' -f3"), emit: versions_hifitrimmer, topic: versions
   tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | sed -e "s/samtools //"'), emit: versions_samtools, topic: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ?: ''
   def args2 = task.ext.args2 ?: ''
   def suffix = args.contains('-f') ? "fastq.gz"  : "fasta.gz"
   def input_convert = !input.name.endsWith('bam') ? "<(samtools import ${input} ${args2} -@ $task.cpus)" : input
   """
   hifi_trimmer filter_bam \\
      -t ${task.cpus} \\
      ${args} \\
      ${input_convert} \\
      ${bed} \\
      ${prefix}.${suffix}
   """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def suffix = args.contains('-f') ? "fastq.gz"  : "fasta.gz"
   """
   echo "stub" | gzip > ${prefix}.${suffix}
   """
}
