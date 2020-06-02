#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Include NfUtils
params.classpath = "umi_tools/groovy/NfUtils.groovy"
Class groovyClass = new GroovyClassLoader(getClass().getClassLoader()).parseClass(new File(params.classpath));
GroovyObject nfUtils = (GroovyObject) groovyClass.newInstance();

// Define internal params
module_name = 'dedup'

// Local default params
params.internal_outdir = 'results'
params.internal_process_name = 'dedup'

// Check for internal parameter overrides
nfUtils.check_internal_overrides(module_name, params)

/*-------------------------------------------------> DEDUP OPTIONS <-----------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------
OUTPUT STATS OPTION -> --output-stats=[PREFIX]
-------------------------------------------------------------------------------------------------------------------------------*/
//Enter the prefix 
params.internal_output_stats = 'sample'

/*-----------------------------------------------------------------------------------------------------------------------------
EXTRACT BARCODE OPTION -> --extract-umi-method
-------------------------------------------------------------------------------------------------------------------------------*/

//Default extract barcode -> Barcodes are contained at the end of the read separated as specified with --umi-separator option
//Enter the character that corresponds to the UMI separator
params.internal_umi_separator = ':'

/*-----------------------------------------------------------------------------------------------------------------------------*/

//Barcodes contained in a tag(s) -> --extract-umi-method=tag
params.internal_extract_method_tag = false

// IF params.internal_extract_method_tag = true -> choose between the 6 options below
//Umi tag options
//--umi-tag=[TAG], --umi-tag-split=[SPLIT], --umi-tag-delimiter=[DELIMITER]
//Select tag, split ot delimiter
params.internal_umi_tag_active = false
params.internal_umi_tag_split_active = false
params.internal_umi_tag_delim_active = false
//If one set to true, insert tag, split or delimiter in the param below
params.internal_umi_tag = ''


//Cell tag options
//--umi-tag=[TAG], --umi-tag-split=[SPLIT], --umi-tag-delimiter=[DELIMITER]
//Select tag, split ot delimiter
params.internal_cell_tag_active = false
params.internal_cell_tag_split_active = false
params.internal_cell_tag_delim_active = false
//If one set to true, insert tag, split or delimiter in the param below
params.internal_cell_tag = ''

/*-----------------------------------------------------------------------------------------------------------------------------
OUTPUT OPTION -> --stdout or -S
-------------------------------------------------------------------------------------------------------------------------------*/
//Insert output file name
params.internal_output_file_name = "sample.dedup.bam"

// dedup reusable component
process dedup {
    publishDir "umi_tools/dedup/${params.internal_outdir}/${params.internal_process_name}",
        mode: "copy", overwrite: true

    input:
      tuple val(sample_id), path(bai), path(bam)
       
    output:
      tuple val(sample_id), path(bam), emit: dedupBam

    script:

    //Initializing dedup arguments
    dedup_pre_args = "umi_tools dedup "
    dedup_post_args = ''

    //Method used to extract barcodes
    if (params.internal_umi_separator != null){
      dedup_pre_args += "--umi-separator=\"$params.internal_umi_separator\" "
    }

    //Input args 
    dedup_pre_args += "-I "

    //Output_stats option
    if (params.internal_output_stats != null){
      dedup_post_args += "--output-stats=$params.internal_output_stats "
    }

    //Output/ stdout option 
    if (params.internal_output_file_name != null){
        dedup_post_args += "-S $params.internal_output_file_name "
    }

     // Displays the umi_tools command line to check for mistakes
    println dedup_pre_args
    println dedup_post_args

    """
    $dedup_pre_args $bam $dedup_post_args
    """
}


//fileName=`basename $bam`
//sampleName="\${fileName%.Aligned.sortedByCoord.out.bam}"
//umi_tools dedup --umi-separator=":" -I $bam -S \${sampleName}.dedup.bam --output-stats=\${sampleName}