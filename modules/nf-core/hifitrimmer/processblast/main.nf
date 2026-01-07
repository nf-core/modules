process HIFITRIMMER_PROCESSBLAST {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/718820d53cc02c7244d2be27e4fe89cd325fff6a3bff2fec0d2d2e4030e7efcf/data' :
      'community.wave.seqera.io/library/hifi_trimmer:2.1.0--7bdb23c108803277' }"

   input:
   tuple val(meta), path(bam), path(blast)
   path(yaml)


   output:
   tuple val(meta), path("*.bed.gz")      , emit: bed
   tuple val(meta), path("*.summary.json"), emit: summary
   tuple val(meta), path("*.hits")        , emit: hits, optional: true
   tuple val("${task.process}"), val('hifi_trimmer'), eval("hifi_trimmer --version | sed 's/hifi_trimmer, version*//'"), emit: versions_hifitrimmer, topic: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ? task.ext.args : ''
   def args1 = task.ext.args1 ? task.ext.args1 : ''
   """
   hifi_trimmer process_blast $args $blast $yaml -t ${task.cpus}
   """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   echo "stub" | gzip > ${prefix}.bed.gz
   echo "stub" | gzip > ${prefix}.summary.json
   echo "stub" | gzip > ${prefix}.hits
   """
}
