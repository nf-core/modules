#!/usr/bin/env nextflow

// Include NfUtils
params.classpath = "umi_tools/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'dedup'

// Specify DSL2
nextflow.preview.dsl = 2

// Local default params
params.internal_outdir = 'results'
params.internal_process_name = 'dedup'

// Check for internal parameter overrides
nfUtils.check_internal_overrides(module_name, params)

// dedup reusable component
process dedup {
    publishDir "umi_tools/dedup/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(bai), path(bam)
       
    output:
      tuple val(sample_id), path(bam), emit: dedupBam

    script:
    """
    fileName=`basename $bam`
    sampleName="\${fileName%.Aligned.sortedByCoord.out.bam}"
    umi_tools dedup --umi-separator=":" -I $bam -S \${sampleName}.dedup.bam --output-stats=\${sampleName}
    """
}