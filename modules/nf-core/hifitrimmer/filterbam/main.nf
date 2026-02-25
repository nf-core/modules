process HIFITRIMMER_FILTERBAM {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2f9ae9ae67cf4c1ee39387e2607cac90cd699a6233ed742970dfc29e5ab43539/data' :
      'community.wave.seqera.io/library/hifi_trimmer_htslib_samtools:db636d1a487afe47' }"

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
   def input_convert = !input.name.endsWith('bam') ? "<(samtools import ${input} ${args2} -@ ${task.cpus})" : input
   """
   hifi_trimmer filter_bam \\
      -t ${task.cpus} \\
      ${args} \\
      ${input_convert} \\
      ${bed} \\
      ${prefix}.${suffix}
   """

   stub:
   def args = task.ext.args ?: ''
   def prefix = task.ext.prefix ?: "${meta.id}"
   def suffix = args.contains('-f') ? "fastq.gz"  : "fasta.gz"
   """
   echo "stub" | gzip > ${prefix}.${suffix}
   echo ${args}
   """
}
