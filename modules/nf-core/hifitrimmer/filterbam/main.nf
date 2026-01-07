process HIFITRIMMER_FILTERBAM {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/718820d53cc02c7244d2be27e4fe89cd325fff6a3bff2fec0d2d2e4030e7efcf/data' :
      'community.wave.seqera.io/library/hifi_trimmer:2.1.0--7bdb23c108803277' }"

   input:
   tuple val(meta), path(bam)
   tuple val(meta1), path(bed)


   output:
   tuple val(meta), path("*.hifi_trimmer.fast{q,a}*"), emit: filtered
   path  "versions.yml"                              , emit: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ? task.ext.args : ''
   def suffix = args.contains('-f') ? "fastq.gz"  : "fasta"
   """
   hifi_trimmer filter_bam $args $bam $bed ${prefix}.hifi_trimmer.${suffix} -t ${task.cpus}

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
      \$(hifi_trimmer --version | sed 's/, version/: /')
   END_VERSIONS
   """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   echo "stub" | gzip > ${prefix}.hifi_trimmer.fastq.gz

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
      \$(hifi_trimmer --version | sed 's/, version/: /')
   END_VERSIONS
   """
}
