process BCL2FASTQ {
    tag "${ meta.lane ? meta.id + "." + meta.lane : meta.id }"
    label 'process_high'

    container "nf-core/bcl2fastq:2.20.0.422"

    input:
    tuple val(meta), path(samplesheet), path(run_dir)

    output:
    tuple val(meta), path("output/**_S[1-9]*_R?_00?.fastq.gz"), emit: fastq
    tuple val(meta), path("output/**_S[1-9]*_I?_00?.fastq.gz"), optional: true, emit: fastq_idx
    tuple val(meta), path("output/**Undetermined_S0*_R?_00?.fastq.gz"), optional: true, emit: undetermined
    tuple val(meta), path("output/**Undetermined_S0*_I?_00?.fastq.gz"), optional: true, emit: undetermined_idx
    tuple val(meta), path("output/Reports"), emit: reports
    tuple val(meta), path("output/Stats"), emit: stats
    tuple val(meta), path("InterOp/*.bin"), emit: interop
    tuple val("${task.process}"), val('bcl2fastq'), eval("bcl2fastq -V 2>&1 | grep -m 1 bcl2fastq | sed 's/^.*bcl2fastq v//'"), topic: versions, emit: versions_bcl2fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("BCL2FASTQ module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def input_tar = run_dir.toString().endsWith(".tar.gz") ? true : false
    def input_dir = input_tar ? run_dir.toString() - '.tar.gz' : run_dir
    """
    if [ ! -d ${input_dir} ]; then
        mkdir -p ${input_dir}
    fi

    if ${input_tar}; then
        ## Ensures --strip-components only applied when top level of tar contents is a directory
        ## If just files or multiple directories, place all in ${input_dir}

        if [[ \$(tar -taf ${run_dir} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
            tar \\
                -C ${input_dir} --strip-components 1 \\
                -xavf \\
                ${args2} \\
                ${run_dir} \\
                ${args3}
        else
            tar \\
                -C ${input_dir} \\
                -xavf \\
                ${args2} \\
                ${run_dir} \\
                ${args3}
        fi
    fi

    bcl2fastq \\
        ${args} \\
        --output-dir output \\
        --runfolder-dir ${input_dir} \\
        --sample-sheet ${samplesheet} \\
        --processing-threads ${task.cpus}

    cp -r ${input_dir}/InterOp .
    """

    stub:
    """
    mkdir -p output
    echo "" | gzip > output/Sample1_S1_L001_R1_001.fastq.gz
    echo "" | gzip > output/Undetermined_S0_L001_R1_001.fastq.gz

    mkdir -p InterOp
    echo "fake interop file" > InterOp/IndexMetricsOut.bin
    mkdir -p output/Reports
    echo "fake report" > output/Reports/index.html
    mkdir -p output/Stats
    echo "fake stats" > output/Stats/Stats.json
    """
}

def generateReadgroup(ch_fastq) {
    ch_fastq
        .transpose()
        .map { fc_meta, fastq ->
            def meta = [:]
            meta.id = fastq.getSimpleName().toString() - ~/_R[0-9]_001.*$/
            meta.samplename = fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/
            meta.fcid = fc_meta.id
            meta.lane = fc_meta.lane
            // The buffered input stream allows reading directly from cloud storage
            // It will not make a local copy of the file.
            def line = ""
            fastq.withInputStream { fq ->
                def gzipStream = new java.util.zip.GZIPInputStream(fq)
                def decoder = new InputStreamReader(gzipStream, 'ASCII')
                def buffered = new BufferedReader(decoder)
                line = buffered.readLine()
                buffered.close()
            }
            if (line != null && line.startsWith('@')) {
                line = line.substring(1)
                // expected format is like:
                // xx:yy:FLOWCELLID:LANE:... (seven fields)
                def fields = line.split(':')
                // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
                // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
                //def sequencer_serial = fields[0]
                //def run_nubmer       = fields[1]
                def fcid = fields[2]
                def lane = fields[3]
                def index = fields[-1] =~ /[GATC+-]/ ? fields[-1] : ""
                def ID = [index, lane].join(".")
                def LB = ""
                def PL = "ILLUMINA"
                def PU = [fcid, lane].findAll().join(".")
                def SM = fastq.getSimpleName().toString() - ~/_S[0-9]+.*$/
                meta.readgroup = ["ID": ID, "SM": SM, "PL": PL, "PU": PU, "LB": LB]
            }
            else {
                println("No reads were found in FASTQ file: ${fastq}")
                meta.readgroup = [:]
            }
            return [meta, fastq]
        }
        .groupTuple(by: [0])
        .map { meta, fastq ->
            meta.single_end = fastq.size() == 1
            return [meta, fastq.flatten()]
        }
}
