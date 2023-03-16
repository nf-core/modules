process CREATESNVS {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each preprocess_data
        val (params)

    output:
        tuple val(vcf), val(rnaID), val(bam_file),
            path("*.vcf.gz"), path("*.vcf.gz.tbi")       , emit: result

    shell:
        vcf     = preprocess_data.vcf
        rnaID   = preprocess_data.rnaID

        vcf_file    = preprocess_data.vcf_file
        bam_file    = preprocess_data.bam_file

        bcftoolsCmd = params.bcftoolsCmd
        samtoolsCmd = params.samtoolsCmd

        ncbi2ucsc   = params.ncbi2ucsc
        ucsc2ncbi   = params.ucsc2ncbi
        '''
            #!/bin/bash
            set -e

            # 1 {input.ncbi2ucsc}
            # 2 {input.ucsc2ncbi}
            # 3 {input.vcf_file}
            # 4 {wildcards.vcf}
            # 5 {input.bam_file}
            # 6 {output.snvs_filename}
            # 7 {config[tools][bcftoolsCmd]}
            # 8 {config[tools][samtoolsCmd]}

            ncbi2ucsc=!{ncbi2ucsc}
            ucsc2ncbi=!{ucsc2ncbi}
            vcf_file=!{vcf_file}
            vcf_id=!{vcf}
            bam_file=!{bam_file}
            output="!{vcf}--!{rnaID}.vcf.gz"
            bcftools=!{bcftoolsCmd}
            samtools=!{samtoolsCmd}

            tmp=$(mktemp)
            tmp2=$(mktemp)

            canonical_chr="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"

            echo 'Filter SNVs'

            # if not doing QC, check for sampleID and select heterozygous variants for MAE
            if [ $vcf_id != 'QC' ]; then
                # match the sampleID from the vcf file
                sample_flag="-s ${vcf_id}"
                sample_name="-sn ${vcf_id}"
                # JEXL pattern to find the heterozygous genotypes
                select_pattern="-select vc.getGenotype('${vcf_id}').isHet()"

            else
                # when doing QC we don't have a match for the sample
                sample_flag=""
                # when doing QC we want all of our QC variants. so don't filter based on GT
                # empty sample_name and select_pattern will not filter the vcf file
                sample_name=""
                select_pattern=""
            fi

            # view the vcf file and remove the info header information and the set the INFO column to '.'
            # split any multi-allelic lines
            # pull out the sample and only the snps that have at least 2 reads supporting it
            # Using bcftools -Ou to speed up processing
            $bcftools view  $vcf_file -r $canonical_chr | \
                grep -vP '^##INFO=' | \
                awk -F'\t' 'BEGIN {OFS = FS} { if($1 ~ /^[^#]/){ $8 = "." }; print $0 }' | \
                $bcftools norm -Ou -m-both | \
                $bcftools norm -Ou -d both | \
                $bcftools view ${sample_flag} -m2 -M2 -v snps > ${tmp2}

            # use the select_pattern defined above to pull out the heterozygous variants used for MAE
            gatk SelectVariants -V ${tmp2} ${sample_name} ${select_pattern} -O ${tmp}

            # explicitly remove all BIALLELIC variants. This is needed as ASEReadCounter does not support them.
            gatk SelectVariants --restrict-alleles-to BIALLELIC -V ${tmp} -O ${tmp2}

            # zip and save as tmp file
            bgzip -c ${tmp2} > ${tmp}
            $bcftools index -t ${tmp}

            rm -f $tmp2

            # compare and correct chromosome format mismatch
            bam_chr=$($samtools idxstats ${bam_file} | cut -f1 | grep "^chr" | wc -l)
            vcf_chr=$($bcftools index --stats $tmp   | cut -f1 | grep "^chr" | wc -l)

            if [ ${vcf_chr} -eq 0  ] && [ ${bam_chr} -ne 0 ]  # VCF: UCSC, BAM: NCBI
            then
                echo "converting from NCBI to UCSC format"
                $bcftools annotate --rename-chrs $ncbi2ucsc $tmp | bgzip > ${output}
                rm ${tmp}
                rm ${tmp}.tbi
            elif [ ${vcf_chr} -ne 0  ] && [ ${bam_chr} -eq 0 ]  # VCF: NCBI, BAM: UCSC
            then
                echo "converting from UCSC to NCBI format"
                $bcftools annotate --rename-chrs $ucsc2ncbi $tmp | bgzip > ${output}
                rm ${tmp}
                rm ${tmp}.tbi
            else  # VCF and BAM have same chromosome format
                mv $tmp ${output}
                rm ${tmp}.tbi
            fi

            num_out=$(zcat "${output}" | grep -vc '#' )
            if [ "${num_out}" -eq 0 ]
            then
            printf  "%s\n" "" "ERROR: No entries after filtering for SNVs" \
            "  Make sure that the VCF file is correctly formatted and contains heterozygous variants." \
            "  This analysis is independent per sample, so consider removing the sample from your analysis as a last resort." \
            "" "  VCF ID: ${vcf_id}" \
            "  VCF file: ${vcf_file}" \
            "  BAM file: ${bam_file}"
            exit 1
            fi

            $bcftools index -t ${output}
        '''
}
