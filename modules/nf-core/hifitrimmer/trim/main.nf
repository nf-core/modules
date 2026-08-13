process HIFITRIMMER_TRIM {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/21/2153598b48eda939e4ac790e5f6e0d28d714faa2293dc04539a9a3f51b5a943d/data' :
      'community.wave.seqera.io/library/hifi_trimmer:5.0.0--c5b9bef0c5bd186f' }"

   input:
   tuple val(meta), path(input), path(bed)


   output:
   tuple val(meta), path("${prefix}.fasta"), emit: fasta, optional: true
   tuple val(meta), path("${prefix}.fastq"), emit: fastq, optional: true
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
   def suffix = args.contains('-f cram') ? 'cram'
      : args.contains('-f sam') ? 'sam'
      : args.contains('-f bam') ? 'bam'
      : args.contains('-f fastq') ? 'fastq' : 'fasta'

   // Check compatibility of input-output format combinations before the process runs
   def inputName = input.name.toLowerCase()
   if (inputName =~ /\.(fa|fasta)(\.gz)?$/ && suffix != 'fasta') {
      error "ERROR: FASTA input can only produce FASTA output, but '-f ${suffix}' was requested."
   }
   if (inputName =~ /\.(fq|fastq)(\.gz)?$/ && suffix in ['bam', 'sam', 'cram']) {
      error "ERROR: FASTQ input cannot produce ${suffix.toUpperCase()} output."
   }
   if (input.name == "${prefix}.${suffix}") {
      error "ERROR: Output file '${prefix}.${suffix}' collides with input file name. Please set a different prefix via task.ext.prefix or meta.id."
   }

   """
   hifi-trimmer trim \\
      -t ${task.cpus} \\
      ${args} \\
      ${input} \\
      ${bed} \\
      > ${prefix}.${suffix}
   """

   stub:
   prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ?: ''
   def suffix = args.contains('-f cram') ? 'cram'
      : args.contains('-f sam') ? 'sam'
      : args.contains('-f bam') ? 'bam'
      : args.contains('-f fastq') ? 'fastq' : 'fasta'
   """
   printf "stub\n" > ${prefix}.${suffix}
   echo ${args}
   """
}
