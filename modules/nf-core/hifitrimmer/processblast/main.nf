process HIFITRIMMER_PROCESSBLAST {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1e/1e5760f3cc4b6cc405353f1122994fd4ca6defd931985f6d0cba3c6ca72e43ab/data' :
      'community.wave.seqera.io/library/hifi_trimmer:2.2.0--1b370153702e2fcc' }"

   input:
   tuple val(meta), path(blast)
   tuple val(meta1), path(yaml)


   output:
   tuple val(meta), path("*.bed.gz")      , emit: bed
   tuple val(meta), path("*.summary.json"), emit: summary
   tuple val(meta), path("*.hits")        , emit: hits, optional: true
   tuple val("${task.process}"), val('hifi_trimmer'), eval("hifi_trimmer --version | cut -d' ' -f3"), emit: versions_hifitrimmer, topic: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ? task.ext.args : ''
   """
   hifi_trimmer process_blast \\
      -t ${task.cpus} \\
      --prefix ${prefix} \\
      ${args} \\
      ${blast} \\
      ${yaml}
   """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   echo "stub" | gzip > ${prefix}.bed.gz
   echo "stub" | gzip > ${prefix}.summary.json
   echo "stub" | gzip > ${prefix}.hits
   """
}
