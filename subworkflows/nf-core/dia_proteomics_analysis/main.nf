/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { DIANN_INSILICOLIBRARYGENERATION } from '../../../modules/nf-core/diann/insilicolibrarygeneration/main'
//include { DIANN_PRELIMINARYANALYSIS      } from '../../../modules/nf-core/diann/preliminaryanalysis/main'
//include { DIANN_ASSEMBLEEMPIRICALLIBRARY } from '../../../modules/nf-core/diann/assembleempiricallibrary/main'
//include { DIANN_INDIVIDUALANALYSIS       } from '../../../modules/nf-core/diann/individualanalysis/main'
//include { DIANN_FINALQUANTIFICATION      } from '../../../modules/nf-core/diann/finalquantification/main'
include { QUANTMSUTILS_DIANNCFG          } from '../../../modules/nf-core/quantmsutils/dianncfg/main'
include { QUANTMSUTILS_MZMLSTATISTICS    } from '../../../modules/nf-core/quantmsutils/mzmlstatistics/main'
include { QUANTMSUTILS_DIANN2MZTAB       } from '../../../modules/nf-core/quantmsutils/diann2mztab/main'
include { MSSTATS_MSSTATSLFQ             } from '../../../modules/nf-core/msstats/msstatslfq/main'

include { DIANN as DIANN_INSILICOLIBRARYGENERATION } from '../../../modules/nf-core/diann/main'
include { DIANN as DIANN_PRELIMINARYANALYSIS } from '../../../modules/nf-core/diann/main'
include { DIANN as DIANN_ASSEMBLEEMPIRICALLIBRARY } from '../../../modules/nf-core/diann/main'
include { DIANN as DIANN_INDIVIDUALANALYSIS } from '../../../modules/nf-core/diann/main'
include { DIANN as DIANN_FINALQUANTIFICATION } from '../../../modules/nf-core/diann/main'

def sortByFilename = { a, b -> file(a).getName() <=> file(b).getName() }

def sortListsByPathName = { tuple, sortIndex ->
    def meta = tuple[0]
    def sortOrder = tuple[sortIndex].withIndex()
        .sort { it[0].name }
        .collect { it[1] }
    
    [meta] + (1..<tuple.size()).collect { i ->
        sortOrder.collect { j -> tuple[i][j] }
    }
}

def extractDiannMassAccuracyFromLog = { diann_log ->
    def settingsLine = diann_log.text.split('\n').find { it.contains('Averaged recommended settings') }
    if (settingsLine) {
        def fields = settingsLine.split()
        return [
            mass_acc_ms1: fields[14]?.replaceAll(/[^0-9]/, ''),
            mass_acc_ms2: fields[10]?.replaceAll(/[^0-9]/, ''),
            scan_window: fields[18]?.replaceAll(/[^0-9]/, '')
        ]
    }
    return [mass_acc_ms1: null, mass_acc_ms2: null, scan_window: null]
}

def defineMassAccuracySettings = { logSettings, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, wf_scan_window, wf_mass_acc_automatic, wf_scan_window_automatic, wf_pg_level ->
    def mass_acc_ms1, mass_acc_ms2, scan_window
    
    if (wf_mass_acc_automatic || wf_scan_window_automatic) {
        mass_acc_ms1 = logSettings.mass_acc_ms1
        mass_acc_ms2 = logSettings.mass_acc_ms2
        scan_window = logSettings.scan_window
    } else if (precursor_tolerance_unit?.toLowerCase().endsWith('ppm') && fragment_tolerance_unit?.toLowerCase().endsWith('ppm')) {
        mass_acc_ms1 = precursor_tolerance
        mass_acc_ms2 = fragment_tolerance
        scan_window = wf_scan_window
    } else {
        mass_acc_ms1 = logSettings.mass_acc_ms1
        mass_acc_ms2 = logSettings.mass_acc_ms2
        scan_window = logSettings.scan_window
    }
    
    [mass_acc_ms1: mass_acc_ms1, mass_acc_ms2: mass_acc_ms2, scan_window: scan_window, pg_level: wf_pg_level]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DIA_PROTEOMICS_ANALYSIS {
    take:
    ch_input                 // Channel of tuples of val(meta), path(ms_file), val(enzyme), val(fixed_mods), val(variable_mods), val(precursor_tolerance), val(fragment_tolerance), val(precursor_tolerance_unit), val(fragment_tolerance_unit)
    ch_searchdb              // Channel of tuples of [val(meta), path(fasta)]
    ch_expdesign             // Channel of tuples of [val(meta), path(expdesign)]
    random_preanalysis       // Tuple of [boolean, integer, integer] for random preanalysis(?), n, seed
    wf_scan_window           // Integer: Scan window for DIA-NN individual analysis
    wf_mass_acc_automatic    // Boolean: Whether to use automatic mass accuracy from preliminary analysis
    wf_scan_window_automatic // Boolean: Whether to use automatic scan window from preliminary analysis
    wf_diann_debug           // Boolean: Enable DIA-NN debug output
    wf_pg_level              // String: Protein group level for DIA-NN analysis
    ch_speclib               // Channel of tuples of val(meta), path(speclib)
    ch_empirical             // tuple val(meta), path(assembly_log), path(empirical_library)

    main:

    ch_versions = Channel.empty()
    (random_preanalysis, random_preanalysis_n, random_preanalysis_seed) = random_preanalysis ?: [false, null, null]

    input = ch_expdesign
        .combine(ch_searchdb)
        .combine(ch_input)
        .multiMap{ 
            meta_exp, 
            expdesign,
            meta_fasta, 
            fasta, 
            meta_input, 
            ms_file, 
            enzyme, 
            fixed_mods, 
            variable_mods, 
            precursor_tolerance, 
            fragment_tolerance, 
            precursor_tolerance_unit, 
            fragment_tolerance_unit ->

            def dia_params = [fragment_tolerance, fragment_tolerance_unit, precursor_tolerance, precursor_tolerance_unit, enzyme, fixed_mods, variable_mods].join(';')
            
            meta_input = meta_input + [
                'enzyme': enzyme,
                'fixed_mods': fixed_mods,
                'variable_mods': variable_mods,
                'precursor_tolerance': precursor_tolerance,
                'fragment_tolerance': fragment_tolerance,
                'precursor_tolerance_unit': precursor_tolerance_unit,
                'fragment_tolerance_unit': fragment_tolerance_unit,
                'dia_params': dia_params
            ]

            def meta_exp_searchdb = meta_exp + meta_fasta + [id: meta_exp.id + '_' + meta_fasta.id]
            def id_enzyme_mods = (enzyme + '_' + fixed_mods + '_' + variable_mods).replaceAll(/[^a-zA-Z0-9_]/, '_')
            def meta_enzyme_mods = meta_input.subMap(['enzyme', 'fixed_mods', 'variable_mods']) + [id: id_enzyme_mods]
            def meta_fasta_enzyme = [id: meta_fasta.id + '_' + meta_enzyme_mods.id]
            def meta_input_fasta = meta_input + meta_fasta + [id: meta_input.id + '_' + meta_fasta.id] + [experiment: meta_exp_searchdb]
            meta_input.experiment = meta_exp_searchdb

            expdesign: [meta_exp_searchdb, expdesign]

            ms_file: [meta_input, meta_input_fasta, ms_file]
            ms_file_fasta: [meta_fasta_enzyme, meta_input_fasta, ms_file]
            search_db_by_exp: [meta_exp_searchdb, fasta]
            search_db_by_enzyme: [meta_enzyme_mods, meta_fasta_enzyme, meta_fasta, fasta]
            enzyme_mods: [meta_enzyme_mods, enzyme, fixed_mods, variable_mods]
            mass_tolerance_settings: [meta_input_fasta, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, ms_file, fasta]
        }

    //
    // MODULE: Generate DIA-NN configuration. Needs to run once each unique enzyme, fixed modifications, and variable modifications combination.
    //
    // Only run config generation and speclib generation if no speclib is provided

    ch_config_input = input.enzyme_mods.unique() // [meta_enzyme_mods, enzyme, fixed_mods, variable_mods]
        .combine(ch_speclib.ifEmpty(null))
        .filter{tuple -> tuple.any { it == null }} 
        .map{ meta, enzyme, fixed_mods, variable_mods, _ ->
            [meta, enzyme, fixed_mods, variable_mods]
        }
    
    QUANTMSUTILS_DIANNCFG(ch_config_input)
    ch_versions = ch_versions.mix(QUANTMSUTILS_DIANNCFG.out.versions)

    //
    // MODULE: In-silico library generation. Needs to run once for every unique config/ FASTA dataset combination.
    //

    ch_insilico_library_input = input.search_db_by_enzyme    // [meta_enzyme_mods, meta_fasta_enzyme, meta_fasta, fasta]
        .combine(QUANTMSUTILS_DIANNCFG.out.diann_cfg, by: 0) // [meta_enzyme_mods, meta_fasta_enzyme, meta_fasta, fasta, cfg_file]
        .unique() // Multiple inputs might have the same config
        .map { meta_enzyme_mods, meta_fasta_enzyme, meta_fasta, fasta, cfg_file -> 
            [meta_fasta_enzyme + [config: cfg_file.text], [], [], fasta, [], []] // Added ms_file_names as []
        }

    DIANN_INSILICOLIBRARYGENERATION(ch_insilico_library_input)
    ch_speclib = ch_speclib.mix(DIANN_INSILICOLIBRARYGENERATION.out.predict_speclib)
    ch_versions = ch_versions.mix(DIANN_INSILICOLIBRARYGENERATION.out.versions)

    // In-silico libraries have been generated for combinations of FASTA and configuration, which 
    //may have been the same over multiple inputs. Use a combine to annotate inputs with in silico libraries.
    ch_fasta_input_with_speclib = ch_speclib
        .map{ meta, speclib -> [meta.findAll { key, value -> key != 'config' }, speclib] } // Filter out config from meta to allow the combine by: 0
        .combine(input.ms_file_fasta, by: 0)
        .map { meta_fasta_enzyme, speclib, meta_input_fasta, ms_file ->
            [meta_input_fasta, ms_file, speclib]
        }
    
    //
    // MODULE: Preliminary analysis
    //

    // We only need to do preliminary analysis and empirical library assembly if no empirical library is provided

    preliminary_branches = ch_fasta_input_with_speclib  
        .join(ch_empirical, remainder: true)
        .filter { tuple ->
            tuple.any { it == null } // ch_empirical was an empty channel
        }
        .map { meta_input_fasta, ms_file, speclib, _ ->
            [meta_input_fasta, ms_file, speclib]
        }
        .branch { 
            random_preanalysis: random_preanalysis
            no_random_preanalysis: no_random_preanalysis
        }

    ch_preliminary_analysis_input = preliminary_branches.no_random_preanalysis
        .mix(
            preliminary_branches.random_preanalysis
                .toSortedList{ a, b -> file(a[1]).getName() <=> file(b[1]).getName() }
                .flatMap()
                .randomSample(random_preanalysis_n, random_preanalysis_seed)
        )

    DIANN_PRELIMINARYANALYSIS(
        ch_preliminary_analysis_input.map { meta, ms_file, speclib ->
            [meta, ms_file, [], [], speclib, []]  // DIANN module input: [meta, ms_files, ms_file_names, fasta, library, quant]
        }
    )
    ch_versions = ch_versions.mix(DIANN_PRELIMINARYANALYSIS.out.versions)

    //
    // MODULE: Assemble empirical library with all inputs from the same experiment + search DB
    //
    ch_empirical_input = ch_preliminary_analysis_input    // [meta_input_fasta, ms_file, speclib]
        .join(DIANN_PRELIMINARYANALYSIS.out.diann_quant)  // [meta_input_fasta, ms_file, speclib, diann_quant]
        .map { [it[0].experiment] + it.drop(1)}           // [meta_exp_searchdb, ms_file, speclib, diann_quant]
        .groupTuple()                                     // [meta_exp_searchdb, [ms_file], [speclib], [diann_quant]]
        .map{ sortListsByPathName(it, 1) }                // [meta_exp_searchdb, [sorted_ms_file], [sorted_speclib], [sorted_diann_quant]]
        .map { meta, ms_files, speclib, diann_quant ->
            [meta, ms_files, [], [], speclib, diann_quant]    // DIANN module input: add ms_file_names and fasta placeholders
        }

    DIANN_ASSEMBLEEMPIRICALLIBRARY(ch_empirical_input)
    ch_versions = ch_versions.mix(DIANN_ASSEMBLEEMPIRICALLIBRARY.out.versions)

    //
    // Derive suggested settings from DIA-NN log
    //

    ch_empirical_output = DIANN_ASSEMBLEEMPIRICALLIBRARY.out.empirical_library // [meta_exp_searchdb, empirical_library]
        .join(DIANN_ASSEMBLEEMPIRICALLIBRARY.out.log)                          // [meta_exp_searchdb, empirical_library, diann_log]
        .mix(ch_empirical)
        .map{ meta_exp_searchdb, empirical_library, diann_log ->
            def logSettings = extractDiannMassAccuracyFromLog(diann_log)
            [meta_exp_searchdb, logSettings, empirical_library]
        }                                                                      // [meta_exp_searchdb, logSettings, empirical_library]

    //
    // Generate mass accuracy settings combining the settings from the log and the settings from the input
    //

    ch_individual_analysis_input = input.mass_tolerance_settings // [meta_input_fasta, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, ms_file, fasta]
        .map{ [it[0].experiment] + it } // [meta_exp_searchdb, meta_input_fasta, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, ms_file, fasta]
        .combine(ch_empirical_output, by: 0) // [meta_exp_searchdb, meta_input_fasta, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, ms_file, fasta, logSettings, empirical_library]
        .map { meta_exp, meta_input_fasta, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, ms_file, fasta, logSettings, empirical_library ->
            def mass_settings = defineMassAccuracySettings(logSettings, precursor_tolerance, fragment_tolerance, precursor_tolerance_unit, fragment_tolerance_unit, wf_scan_window, wf_mass_acc_automatic, wf_scan_window_automatic, wf_pg_level)
            def meta_with_settings = meta_input_fasta + mass_settings
            [meta_with_settings, ms_file, fasta, empirical_library]  
        }
        
    //
    // MODULE: Individual analysis
    //

    DIANN_INDIVIDUALANALYSIS(
        ch_individual_analysis_input.map { meta, ms_file, fasta, empirical_library ->
            [meta, ms_file, [], fasta, empirical_library, []]  // DIANN module input: [meta, ms_files, ms_file_names, fasta, library, quant]
        }
    )
    ch_versions = ch_versions.mix(DIANN_INDIVIDUALANALYSIS.out.versions)

    //
    // MODULE: Final quantification
    //

    ch_final_quantification_input = ch_individual_analysis_input // [meta_with_settings, ms_file, fasta, empirical_library]
        .join(DIANN_INDIVIDUALANALYSIS.out.diann_quant)          // [meta_with_settings, ms_file, fasta, empirical_library, diann_quant]
        .map { [it[0].experiment] + it.drop(1)}                  // [meta_exp_searchdb, ms_file, fasta, empirical_library, diann_quant]
        .groupTuple()                                            // [meta_exp_searchdb, [ms_file], [fasta], [empirical_library], [diann_quant]]
        .map{ sortListsByPathName(it, 1) }                       // [meta_exp_searchdb, [sorted_ms_file], [sorted_fasta], [sorted_empirical_library], [sorted_diann_quant]]
        .map { meta, ms_files, fasta, empirical_library, diann_quant ->
            [meta, [], ms_files.collect{ it.name}, fasta, empirical_library, diann_quant] // Use ms_file_names instead of ms_files for final quant
        }

    DIANN_FINALQUANTIFICATION(ch_final_quantification_input)
    ch_versions = ch_versions.mix(DIANN_FINALQUANTIFICATION.out.versions)
    ch_final_quantification_combined_output = DIANN_FINALQUANTIFICATION.out.main_report
        .mix(DIANN_FINALQUANTIFICATION.out.report_parquet)       // [meta_exp_searchdb, report]
        .last()                                                  // [meta_exp_searchdb, report]
        .join(DIANN_FINALQUANTIFICATION.out.pg_matrix)           // [meta_exp_searchdb, report, pg_matrix]
        .join(DIANN_FINALQUANTIFICATION.out.pr_matrix)           // [meta_exp_searchdb, report, pg_matrix, pr_matrix]
        .combine(DIANN_FINALQUANTIFICATION.out.versions.last())  // [meta_exp_searchdb, report, pg_matrix, pr_matrix, versions]

    //
    // MODULE: Generate mzML statistics
    //
    ch_mzml_stats_input = input.ms_file
        .map {meta_input, meta_input_fasta, ms_file ->
            [meta_input, ms_file]
        }
        .unique()

    QUANTMSUTILS_MZMLSTATISTICS(ch_mzml_stats_input)
    ch_versions = ch_versions.mix(QUANTMSUTILS_MZMLSTATISTICS.out.versions)

    ch_statistics = QUANTMSUTILS_MZMLSTATISTICS.out.ms_statistics // [meta_input, ms_statistics]
        .map{[it[0].experiment, it[1]]}                           // [meta_exp_searchdb, ms_statistics]
        .groupTuple()                                             // [meta_exp_searchdb, [ms_statistics]]

    //
    // MODULE: Convert results
    //

    // We need some DIAN-NN parameters to convert the results to mzTab. This just takes the first 
    // config from the input, which is what quantms does, but we should consider if these need to 
    // be input-wise values

    ch_first_config = input.ms_file
        .map { meta_input, meta_input_fasta, ms_file ->
            [meta_input.experiment, meta_input]
        }
        .first()

    ch_diann2mztab_input = ch_first_config
        .join(ch_final_quantification_combined_output)             // [meta_exp_searchdb, meta_input, report, pg_matrix, pr_matrix, versions]
        .join(input.expdesign.unique())                            // [meta_exp_searchdb, meta_input, report, pg_matrix, pr_matrix, versions, exp_design]
        .join(ch_statistics)                                       // [meta_exp_searchdb, meta_input, report, pg_matrix, pr_matrix, versions, exp_design, ms_statistics]
        .join(input.search_db_by_exp)                              // [meta_exp_searchdb, meta_input, report, pg_matrix, pr_matrix, versions, exp_design, ms_statistics, fasta]
        .map { meta_exp_searchdb, meta_input, report, pg_matrix, pr_matrix, versions, exp_design, ms_statistics, fasta ->
            [meta_input + meta_exp_searchdb, report, pg_matrix, pr_matrix, versions, exp_design, ms_statistics, fasta]
        }

    QUANTMSUTILS_DIANN2MZTAB(ch_diann2mztab_input)
    ch_versions = ch_versions.mix(QUANTMSUTILS_DIANN2MZTAB.out.versions)

    emit:
    versions                = ch_versions
    diann_report            = DIANN_FINALQUANTIFICATION.out.main_report
    diann_report_parquet    = DIANN_FINALQUANTIFICATION.out.report_parquet
    mzml_statistics         = QUANTMSUTILS_MZMLSTATISTICS.out.ms_statistics
    msstats_in              = QUANTMSUTILS_DIANN2MZTAB.out.out_msstats
    out_triqler             = QUANTMSUTILS_DIANN2MZTAB.out.out_triqler
    final_result            = QUANTMSUTILS_DIANN2MZTAB.out.out_mztab
}