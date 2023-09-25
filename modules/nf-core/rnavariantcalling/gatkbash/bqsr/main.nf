process CHANGE_HEADERS {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each change_headers

    output:
        tuple val(sample), path("*_dropHeader.bam"), path("*_dropHeader.bam.bai"), path("*_newDropHeader.txt"), emit: result

    shell:
        sample      = change_headers.sample
        bam_file    = change_headers.BAM
        bai_file    = change_headers.BAI
        '''
            #!/bin/bash

            # helper file for RVC snakemake rule changeHeader

            # 1 {input.bam}
            # 2 {input.bai}
            # 3 {input.sample}
            # 4 {log}
            # 5 {output.bam}
            # 6 {output.bai}
            # 7 {output.newHeader}
            input_bam=!{bam_file}
            input_bai=!{bai_file}
            sample=!{sample}
            log=!{sample}.log.out
            output_bam=!{sample}_dropHeader.bam
            output_bai=!{sample}_dropHeader.bam.bai
            output_newHeader=!{sample}_newDropHeader.txt


            SM_internalHeader=""
            PL_internalHeader=""
            samtools view -H $input_bam |grep "^@RG" |grep SM |head -1|
            while read header ; do
                for i in $header; do
                    if [[ $i == "SM:"* ]]; then
                        SM_internalHeader=${i:3}
                    fi
                    if [[ $i == "PL:"* ]]; then
                        PL_internalHeader=${i:3}
                    fi
                done

                if [[ $SM_internalHeader == $sample && $PL_internalHeader != "" ]]; then
                    echo "Internal Header $SM_internalHeader matches $sample" |tee $log
                    echo "Internal Header is designated: $SM_internalHeader"  |tee -a $log
                    echo "SampleID is $sample" |tee -a $log
                    ln -f  $input_bam $output_bam
                    ln -f  $input_bai $output_bai
                    echo "Done Linking files"
                    samtools view -H $input_bam > $output_newHeader
                else
                    if [[ $SM_internalHeader != $sample ]]; then
                        echo "WARNING"
                        echo "Internal Header is designated: $SM_internalHeader" |tee $log
                        echo "SampleID is $sample"  |tee -a $log
                        echo "Forcing $SM_internalHeader to match $sample" |tee -a $log

                        # sed using regEx in place substitiute 'SM:' followed by any thing that isn't tab or newLine. and then replace it with the sampleID and the delimiter (tab or newLine) that matched in the 1st expression.
                        samtools view -H $input_bam > $output_newHeader
                        sed -E -i "s/(SM:[^\\t|\\n]*)(\t|\\n*)/SM:${sample}\\2/" ${output_newHeader}
                    fi
                    if [[ $PL_internalHeader == "" ]]; then
                        echo "WARNING"
                        echo "Internal PL Header is not designated: $PL_internalHeader" |tee $log
                        echo "PL cannot be empty. Using a Dummy: dummyPL "  |tee -a $log
                        echo "Forcing $PL_internalHeader to dummyPL" |tee -a $log

                        # sed using regEx in place substitiute '@RG' with @RG\tPL:dummyPL
                        samtools view -H $input_bam > $output_newHeader
                        sed -E -i "s/^@RG/@RG\\tPL:dummyPL/" ${output_newHeader}
                    fi
                    samtools reheader $output_newHeader $input_bam > $output_bam
                    samtools index -b $output_bam
                fi
                echo "new header can be found here:$output_newHeader" |tee -a $log
            done
        '''
}
