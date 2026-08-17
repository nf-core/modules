process HIFITRIMMER_TRIM {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/653461311faa0a2175b75d08405548c1a9b8e704c6b6b70b133c4970cb3d7e5c/data' :
      'community.wave.seqera.io/library/hifi_trimmer_gzip:e72290a9b8b5d246' }"

   input:
   tuple val(meta), path(input), path(bed)


   output:
   tuple val(meta), path("${prefix}.fasta.gz"), emit: fasta, optional: true
   tuple val(meta), path("${prefix}.fastq.gz"), emit: fastq, optional: true
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
      : args.contains('-f fastq') ? 'fastq.gz' : 'fasta.gz'

   // Check compatibility of input-output format combinations before the process runs
   def inputName = input.name.toLowerCase()
   if (inputName =~ /\.(fa|fasta)(\.gz)?$/ && suffix != 'fasta.gz') {
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
      -o ${prefix}.${suffix} \\
      ${input} \\
      ${bed}
   """

   stub:
   prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ?: ''
   def suffix = args.contains('-f cram') ? 'cram'
      : args.contains('-f sam') ? 'sam'
      : args.contains('-f bam') ? 'bam'
      : args.contains('-f fastq') ? 'fastq.gz' : 'fasta.gz'
   """
   ${suffix.endsWith('.gz') ? "echo \"\" | gzip > ${prefix}.${suffix}" : "echo \"\" > ${prefix}.${suffix}"}
   echo ${args}
   """
}
