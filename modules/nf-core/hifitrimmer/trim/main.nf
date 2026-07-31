process HIFITRIMMER_TRIM {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/72/720bfb6a0c816137d958bc8658cbfd184f0a41e3e69d43e9fcdb93275620d128/data' :
      'community.wave.seqera.io/library/hifi_trimmer_htslib_samtools:a06b8fcc843c1e68' }"

   input:
   tuple val(meta), path(input), path(bed)


   output:
   tuple val(meta), path("${prefix}.sam"), emit: sam, optional: true
   tuple val(meta), path("${prefix}.bam"), emit: bam, optional: true
   tuple val(meta), path("${prefix}.cram"), emit: cram, optional: true
   tuple val("${task.process}"), val('hifi-trimmer'), eval("hifi-trimmer --version | cut -d' ' -f3"), emit: versions_hifitrimmer, topic: versions
   tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | sed -e "s/samtools //"'), emit: versions_samtools, topic: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ?: ''
   def args2 = task.ext.args2 ?: ''
   def args3 = task.ext.args3 ?: ''
   def suffix = args.contains('-f cram') ? 'cram' : args.contains('-f sam') ? 'sam' : 'bam'
   def input_convert = input.name.endsWith('cram') ? "<(samtools view ${input} -u ${args3} -@ ${task.cpus})" :
        !input.name.endsWith('bam') ? "<(samtools import ${input} ${args2} -@ ${task.cpus})" : input

   if (input.name == "${prefix}.${suffix}") {
      error "ERROR: Output file '${prefix}.${suffix}' collides with input file name. Please set a different prefix via task.ext.prefix or meta.id."
   }

   """
   hifi-trimmer trim \\
      -t ${task.cpus} \\
      ${args} \\
      ${input_convert} \\
      ${bed} \\
      > ${prefix}.${suffix}
   """

   stub:
   prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ?: ''
   def suffix = args.contains('-f cram') ? 'cram' : args.contains('-f sam') ? 'sam' : 'bam'
   """
   printf "stub\n" > ${prefix}.${suffix}
   echo ${args}
   """
}
