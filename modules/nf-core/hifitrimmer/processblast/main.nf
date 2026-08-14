process HIFITRIMMER_PROCESSBLAST {
   tag "$meta.id"
   label 'process_medium'

   conda "${moduleDir}/environment.yml"
   container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/653461311faa0a2175b75d08405548c1a9b8e704c6b6b70b133c4970cb3d7e5c/data' :
      'community.wave.seqera.io/library/hifi_trimmer_gzip:e72290a9b8b5d246' }"

   input:
   tuple val(meta), path(blast)
   tuple val(meta2), path(yaml)


   output:

   tuple val(meta), path("*.bed.gz")      , emit: bed
   tuple val(meta), path("*.summary.json"), emit: summary
   tuple val(meta), path("*.hits")        , emit: hits, optional: true
   tuple val("${task.process}"), val('hifi-trimmer'), eval("hifi-trimmer --version | cut -d' ' -f3"), emit: versions_hifitrimmer, topic: versions

   when:
   task.ext.when == null || task.ext.when

   script:
   def prefix = task.ext.prefix ?: "${meta.id}"
   def args = task.ext.args ? task.ext.args : ''
   """
   hifi-trimmer process-blast \\
      -t ${task.cpus} \\
      --prefix ${prefix} \\
      ${args} \\
      ${blast} \\
      ${yaml}
   """

   stub:
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   echo "" | gzip > ${prefix}.bed.gz
   touch ${prefix}.summary.json
   touch ${prefix}.hits
   """
}
