#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ASHLAR } from '../../../../../modules/modules/nf-core/ashlar/main.nf'
// we zero out the UUID of output tiff images with ZERO_UUID so we get a consistent md5sum
include { ZERO_UUID } from './zero_uuid.nf'

def TEST_ORIG_IMG = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-0{6,7}.ome.tiff"
def TEST_ORIG_IMG_1 = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-06.ome.tiff"
def TEST_ORIG_IMG_2 = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/raw/exemplar-001-cycle-07.ome.tiff"
def TEST_ORIG_DFP_1 = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/illumination/exemplar-001-cycle-06-dfp.tif"
def TEST_ORIG_FFP_1 = "/home/pollen/HITS/nextflow/mcmicro/exemplar-001/illumination/exemplar-001-cycle-06-ffp.tif"

def TEST_IMG_ALL = "/home/pollen/HITS/nextflow/mcmicro/cycif-tonsil/cycif-tonsil-cycle{1,2,3}.ome.tif"
def TEST_IMG_1 = "/home/pollen/HITS/nextflow/mcmicro/cycif-tonsil/cycif-tonsil-cycle1.ome.tif"
def TEST_IMG_2 = "/home/pollen/HITS/nextflow/mcmicro/cycif-tonsil/cycif-tonsil-cycle2.ome.tif"
def TEST_IMG_3 = "/home/pollen/HITS/nextflow/mcmicro/cycif-tonsil/cycif-tonsil-cycle3.ome.tif"

workflow test_ashlar_all_files {

    input_list =  [ [ [ id:'test_all', args: '' ],
               "${TEST_IMG_ALL}" ] ]
    input_channel = Channel.fromList(input_list)

    ASHLAR ( input_channel )

    ZERO_UUID ( ASHLAR.out[0], "25169643" )

}

