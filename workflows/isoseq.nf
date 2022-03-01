/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIsoseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.primers ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// MODULE: Local to the pipeline
//
include { PERL_BIOPERL }        from '../modules/local/perl/bioperl/main'
include { GSTAMA_FILELIST }     from '../modules/local/gstama/filelist/main'
include { GSTAMA_POLYACLEANUP } from '../modules/local/gstama/polyacleanup/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
//include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { PBCCS }            from '../modules/nf-core/modules/pbccs/main'
include { LIMA }             from '../modules/nf-core/modules/lima/main'
include { ISOSEQ3_REFINE }   from '../modules/nf-core/modules/isoseq3/refine/main'
include { BAMTOOLS_CONVERT } from '../modules/nf-core/modules/bamtools/convert/main'
include { MINIMAP2_ALIGN }   from '../modules/nf-core/modules/minimap2/align/main'
include { ULTRA_PIPELINE }   from '../modules/nf-core/modules/ultra/pipeline/main'
include { SAMTOOLS_SORT }    from '../modules/nf-core/modules/samtools/sort/main'
include { GSTAMA_COLLAPSE }  from '../modules/nf-core/modules/gstama/collapse/main'
include { GSTAMA_MERGE }     from '../modules/nf-core/modules/gstama/merge/main'
include { MULTIQC }          from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ISOSEQ {
    //
    // SET UP CHANNELS
    //

    ch_versions = Channel.empty()
    // Check if input directory contains bam and pbis files and set up channels
    if (params.input) {
        Channel // Check presence of bam files
        .fromPath(params.input + '/*.bam')
        .ifEmpty { exit 1, "Cannot find any bam(s) in input directory: ${params.input}\nNB: File names must finish by '.bam'\nNB: Path needs to be enclosed in quotes!" }
        Channel // Check presence of bam.pbi files
        .fromPath(params.input + '/*.bam.pbi')
        .ifEmpty { exit 1, "Cannot find any pbi(s) in input directory: ${params.input}\nNB: File names must finish by '.bam.pbi'\nNB: Path needs to be enclosed in quotes!" }

        Channel // Prepare channel for bams: duplicates item by num of chunk and extract id
        .fromPath(params.input + '/*.bam')
        .flatMap {  // Duplicate each samples to match to the number of chunks
            def array = []
            for ( i = 1 ; i <= params.chunk ; i++ ) { array << it }
            return array
        }
        .map {  // Set channel to [meta, bam, pbi]
            [
                [ id:it.toString().replaceAll(/.*\\//, '').replaceAll(/.bam$/, '') ],
                file(it),
                file(it.toString().replaceAll(/.bam$/, '.bam.pbi'))
            ]
        }
        .set { ch_pbccs_in }
    } else { exit 1, 'OPTION ERROR: bam/bam.pbi directory not provided or cannot be found!' }

    n_samples = new File(params.input).listFiles().count { it.name ==~ /.*.bam$/ }

    Channel // Prepare the pbccs chunk_num channel
        .from((1..params.chunk).step(1).toList()*n_samples)
        .set { ch_chunk_num }


    if (params.primers) {
        Channel // Prepare value channel for primers used for the library preparation
            .value(file(params.primers))
            .set { ch_primers }
    } else { exit 1, 'OPTION ERROR: primers file not provided or cannot be found!' }


    if (params.fasta) {
        Channel // Prepare value channel for reference genome fasta file => minimap2/uLTRA
            .value(file(params.fasta))
            .set { ch_fasta }
    } else { exit 1, 'OPTION ERROR: fasta file not provided or cannot be found' }


    if (params.ultra == true) {
        Channel // --> Prepare gtf value channel for ultra
            .value(file(params.gtf))
            .set { ch_gtf }
    } else { exit 1, "OPTION ERROR: gtf file not provided or cannot be found" }


    //
    // START PIPELINE
    //
    PBCCS(ch_pbccs_in, ch_chunk_num, params.chunk) // Generate CCS from raw reads

    PBCCS.out.bam // Update meta, update id (+chunkX) and store former id
    .map {
        def chk       = (it[1] =~ /.*\.(chunk\d+)\.bam/)[ 0 ][ 1 ]
        def id_former = it[0].id
        def id_new    = it[0].id + "." + chk
        return [ [id:id_new, id_former:id_former], it[1] ]
    }
    .set { ch_pbccs_bam_updated }

    LIMA(ch_pbccs_bam_updated, ch_primers)   // Remove primers from CCS
    ISOSEQ3_REFINE(LIMA.out.bam, ch_primers) // Discard CCS without polyA tails, remove it from the other
    BAMTOOLS_CONVERT(ISOSEQ3_REFINE.out.bam)
    GSTAMA_POLYACLEANUP(BAMTOOLS_CONVERT.out.data)

    // Align CCS (no cluster path) or singletons + transcripts (cluster path)
    // User can choose between minimap2 and uLTRA aligners
    if (params.ultra == true) {
        ULTRA_PIPELINE(GSTAMA_POLYACLEANUP.out.fasta, ch_fasta, ch_gtf)
        PERL_BIOPERL(ULTRA_PIPELINE.out.sam) // Remove remove reads ending with GAP (N) in CIGAR string
    }
    else {
        MINIMAP2_ALIGN(GSTAMA_POLYACLEANUP.out.out, ch_fasta) // Align read against genome
        PERL_BIOPERL(MINIMAP2_ALIGN.out.paf)               // Remove remove reads ending with GAP (N) in CIGAR string
    }

    SAMTOOLS_SORT(PERL_BIOPERL.out.out)   // Sort and convert sam to bam
    GSTAMA_COLLAPSE(SAMTOOLS_SORT.out.bam, ch_fasta) // Clean gene models

    GSTAMA_COLLAPSE.out.bed // replace id with the former sample id and group files by sample
    .map { [ [id:it[0].id_former], it[1] ] }
    .groupTuple()
    .set { ch_tcollapse }

    GSTAMA_FILELIST(ch_tcollapse, Channel.value("capped"), Channel.value("1,1,1")) // Generate the filelist file needed by TAMA merge

    ch_tcollapse // Synchronized bed files produced by TAMA collapse with file list file generated by GSTAMA_FILELIST
    .join( GSTAMA_FILELIST.out.tsv )
    .set { ch_tmerge_in }

    GSTAMA_MERGE(ch_tmerge_in.map { [ it[0], it[1] ] }, ch_tmerge_in.map { it[2] }) // Merge bed files from one sample into one

    //
    // MODULE: Pipeline reporting
    //
    ch_versions = ch_versions.mix(PBCCS.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(LIMA.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(ISOSEQ3_REFINE.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(BAMTOOLS_CONVERT.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(GSTAMA_COLLAPSE.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(PERL_BIOPERL.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(GSTAMA_MERGE.out.versions.first().ifEmpty(null))
    ch_versions = ch_versions.mix(GSTAMA_POLYACLEANUP.out.versions.first().ifEmpty(null))

    if (params.ultra == true) {
        ch_versions = ch_versions.mix(ULTRA_PIPELINE.out.versions.first().ifEmpty(null))
    }
    else {
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first().ifEmpty(null))
    }

    ch_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_versions }
//  //
//  // SUBWORKFLOW: Read in samplesheet, validate and stage input files
//  //
//  INPUT_CHECK (
//      ch_input
//  )
//  ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
//
//  //
//  // MODULE: Run FastQC
//  //
//  FASTQC (
//      INPUT_CHECK.out.reads
//  )
//  ch_versions = ch_versions.mix(FASTQC.out.versions.first())
//
//  CUSTOM_DUMPSOFTWAREVERSIONS (
//      ch_versions.unique().collectFile(name: 'collated_versions.yml')
//  )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowIsoseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
//  ch_multiqc_files = ch_multiqc_files.mix(PBCCS.out.report_json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(LIMA.out.summary.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(LIMA.out.counts.collect{it[1]}.ifEmpty([]))
//  ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//  ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
