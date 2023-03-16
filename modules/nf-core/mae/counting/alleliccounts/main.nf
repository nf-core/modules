process ALLELIC_COUNTS {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each create_snv
        val (params)

    output:
        tuple val(vcf), val(rnaID), path("*.csv.gz")     , emit: result

    shell:
        vcf     = create_snv.vcf
        rnaID   = create_snv.rnaID

        vcf_file    = create_snv.SNV
        bam_file    = create_snv.BAM

        fasta       = create_snv.fasta

        bcftoolsCmd = params.bcftoolsCmd
        samtoolsCmd = params.samtoolsCmd
        gatkCmd     = params.gatkCmd

        ncbi2ucsc   = params.ncbi2ucsc
        ucsc2ncbi   = params.ucsc2ncbi

        gatkIgnoreHeaderCheck = params.gatkIgnoreHeaderCheck
        '''
            #!/bin/bash

            # 1 {input.ncbi2ucsc}
            # 2 {input.ucsc2ncbi}
            # 3 {input.vcf_file}
            # 4 {input.bam_file}
            # 5 {wildcards.vcf}--{wildcards.rna}
            # 6 {input.fasta}
            # 7 {config[mae][gatkIgnoreHeaderCheck]}
            # 8 {output.counted}
            # 9 {params.bcftools}
            #10 {params.samtools}
            #11 {params.gatk}

            ncbi2ucsc=!{ncbi2ucsc}
            ucsc2ncbi=!{ucsc2ncbi}
            vcf_file=!{vcf_file}
            bam_file=!{bam_file}
            mae_id="!{vcf}--!{rnaID}"
            fasta=!{fasta}
            sanity=!{gatkIgnoreHeaderCheck}
            output="!{vcf}--!{rnaID}.csv.gz"
            bcftools=!{bcftoolsCmd}
            samtools=!{samtoolsCmd}
            gatk=!{gatkCmd}

            tmp=$(mktemp)
            header="contig\\tposition\\tvariantID\\trefAllele\\taltAllele\\t"
            header+="refCount\\taltCount\\ttotalCount\\tlowMAPQDepth\\t"
            header+="lowBaseQDepth\\trawDepth\\totherBases\\timproperPairs"
            echo -e $header >> $tmp

            # get chr format
            bam_chr=$($samtools idxstats ${bam_file} | grep -vP "\t0\t0" | cut -f1 | sort -u) # only chr with coverage
            vcf_chr=$($bcftools view ${vcf_file} | cut -f1 | grep -v '#' | sort -u)
            if [ "$(echo ${vcf_chr} | grep -c 'chr')" -eq 0 ]; then
            echo "use NCBI format"
            canonical=$ncbi2ucsc
            else
            echo "use UCSC format"
            canonical=$ucsc2ncbi
            fi

            # subset to standard chromosomes
            chr_subset=$(comm -12 <(cut -f1 -d" " ${canonical} | sort -u) <(echo "${vcf_chr}"))
            chr_subset=$(comm -12 <(echo "${bam_chr}") <(echo "${chr_subset}") | uniq)

            # ASEReadCounter fails without read groups (RG), this snippet checks for RG in bam file
            # and if RG tag isn't present, lets the user know how to fix it
            if samtools view -H ${bam_file} | grep -q "@RG";then
            printf "BAM contains read groups (RG), continuing with ASEReadCounter...\n"
            else
            printf "%s\n" "" "ERROR: BAM file doesn't contain Read Group Tag (RG)" \
            " RG doesn't exist, it can be added using -" \
            "   gatk AddOrReplaceGroups -R /path/to/reference -I /your/input.bam -O /your/output.bam --QUIET true" \
            " https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-" \
            " Try rerunning this module using the BAM with RG tags"
            exit 1
            fi


            for chr in $chr_subset; do
                $gatk ASEReadCounter -R ${fasta} -I ${bam_file} -V ${vcf_file} -L ${chr} --verbosity ERROR --QUIET true --disable-sequence-dictionary-validation ${sanity} | tail -n+2 >> $tmp
            done

            cat $tmp | awk -v id="${mae_id}" -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "ID"} NR>1{print $0, id}' | bgzip >${output}
            rm ${tmp}

            num_out=$(zcat "${output}" | wc -l )
            if [ "${num_out}" -lt 2 ]
            then
            printf  "%s\n" "" "ERROR: No allele-specific counts" \
                "  Make sure that the chromosome styles of the FASTA reference and BAM file match." \
                "  If that isn't the issue, check that your VCF and BAM files are correctly formatted." \
                "  If this problem persists and if this is your only sample causing issues, consider removing it from your analysis, as a last resort." \
                "" "  MAE ID: ${mae_id}" \
                "  VCF file: ${vcf_file}" \
                "  BAM file: ${bam_file}" \
                "  FASTA file: ${fasta}" \
                " Additionally the ReadGroups may be poorly formed. Please refer to https://gagneurlab-drop.readthedocs.io/en/latest/help.html for more information "
            exit 1
            fi

            zcat ${output} | head
        '''
}
