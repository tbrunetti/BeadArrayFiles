from IlluminaBeadArrayFiles import *
import struct
from io import BytesIO
import os
import sys
import argparse
import pandas
import write_gtc


def manipulate_gtc(bpm, gtcDir, snpsToUpdate, outDir):
    def getGtcInfo(gtc):
        genotype_calls = GenotypeCalls(gtc)
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_AUTOCALL_DATE] = genotype_calls.get_autocall_date(
            )  # key:201
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_AUTOCALL_VERSION] = genotype_calls.get_autocall_version(
            )  # key:300
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_B_ALLELE_FREQS] = genotype_calls.get_ballele_freqs(
            )  # key:1012 - per SNP on all SNPs on chip
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_BASE_CALLS] = genotype_calls.get_base_calls(
            )  # key:1003 - per SNP on all SNPs on chip (options: A/C/T/G/-/I/D)
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CALL_RATE] = genotype_calls.get_call_rate(
            )  # key:1006 - per sample -- total number of valid SNPs with genotype calls divided by total SNPs clustered.  If a SNP has been zeroed out in the cluster file, it is not considered in the tatal # clustered snps; confirmed same value as in BCP
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CLUSTER_FILE] = genotype_calls.get_cluster_file(
            )  # key:100
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CONTROLS_X] = genotype_calls.get_control_x_intensities(
            )  # key:500 - 92 values
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CONTROLS_Y] = genotype_calls.get_control_y_intensities(
            )  # key:501 - 92 values
        data[GenotypeCalls._GenotypeCalls__ID_GC10] = genotype_calls.get_gc10(
        )  # key:1009 - per samples - confirmed same as BCP
        data[GenotypeCalls._GenotypeCalls__ID_GC50] = (
            genotype_calls.get_gc50(), genotype_calls.get_num_calls(),
            genotype_calls.get_num_no_calls(),
            genotype_calls.get_num_intensity_only())  # key:1011 - per sample
        data[GenotypeCalls.
             _GenotypeCalls__ID_GENDER] = genotype_calls.get_gender(
             )  # key:1007 - calculated from chip NOT from sample manifest
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_GENOTYPE_SCORES] = genotype_calls.get_genotype_scores(
            )  # key:1004
        data[GenotypeCalls.
             _GenotypeCalls__ID_GENOTYPES] = genotype_calls.get_genotypes(
             )  # key:1002
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_IMAGING_DATE] = genotype_calls.get_imaging_date(
            )  # key:200
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_LOGR_DEV] = genotype_calls.get_logr_dev(
            )  # key:1008 - the standard deviation of the log(r ratio) across all snps. Essentially estimates noise per sample; confirmed same as in BCP
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_NORMALIZATION_TRANSFORMS] = genotype_calls.get_normalization_transforms(
            )  # key:400
        data[GenotypeCalls.
             _GenotypeCalls__ID_NUM_SNPS] = genotype_calls.get_num_snps(
             )  # key:1
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_PERCENTILES_X] = genotype_calls.get_percentiles_x(
            )  # key:1014 - 3 values
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_PERCENTILES_Y] = genotype_calls.get_percentiles_y(
            )  # key:1015 - 3 values
        data[GenotypeCalls.
             _GenotypeCalls__ID_PLOIDY] = genotype_calls.get_ploidy(
             )  # key:2 - per sample
        data[GenotypeCalls.
             _GenotypeCalls__ID_PLOIDY_TYPE] = genotype_calls.get_ploidy_type(
             )  # key:3 - per sample
        data[GenotypeCalls.
             _GenotypeCalls__ID_RAW_X] = genotype_calls.get_raw_x_intensities(
             )  # key:1000 - per SNP
        data[GenotypeCalls.
             _GenotypeCalls__ID_RAW_Y] = genotype_calls.get_raw_y_intensities(
             )  # key:1001 - per SNP
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_NAME] = genotype_calls.get_sample_name(
             )  # key:10
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SAMPLE_PLATE] = genotype_calls.get_sample_plate(
            )  # key:11
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_WELL] = genotype_calls.get_sample_well(
             )  # key:12
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SCANNER_DATA] = genotype_calls.get_scanner_data(
            )  # key:1005
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SLIDE_IDENTIFIER] = genotype_calls.get_slide_identifier(
            )  # key:1016 - Illumina barcode (not the RxCx part)
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SNP_MANIFEST] = genotype_calls.get_snp_manifest(
            )  # key:101 - name of cluster file used
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_LOGR_RATIOS] = genotype_calls.get_logr_ratios(
            )  # key:1013 per SNP, the log(R ratio), 0 = perfect single copy while above and below indicate copy number anomoalies

        return data

    def updateMetaData(data, metaData):
        import itertools

        dataDict = {}
        metaDataUpdates = metaData.rstrip().split(',')
        for update in metaDataUpdates:
            if update.rstrip().split('=')[0] == 'sampleName':
                dataDict[10] = update.rstrip().split('=')[1]
            elif update.rstrip().split('=')[0] == 'sentrixBarcode':
                dataDict[1016] = update.rstrip().split('=')[1]
            elif update.rstrip().split('=')[0] == 'plateName':
                dataDict[11] = update.rstrip().split('=')[1]
            elif update.rstrip().split('=')[0] == 'well':
                dataDict[12] = update.rstrip().split('=')[1]
            else:
                print(
                    'MetaData {} does not exist; please make sure spelling is correct and case sensitive!  Ignoring...'
                    .format(update.rstrip().split('=')[0]))
                sys.stdout.flush()

        for key, value in dataDict.items():
            data[key] = value

        return data

    def snpUpdate(data, line):
        loc = manifest.names.index(line.rstrip().split()[0])
        originalSnp = data[1003][loc]
        data[1003][loc] = str(line.rstrip().split()[1])
        if ((str(line.rstrip().split()[1])[0] != str(
                line.rstrip().split()[1])[1])
                and (str(line.rstrip().split()[1])[0] != '-')):
            data[1002][loc] = 2
        elif (str(line.rstrip().split()[1])[0] == '-') and (str(
                line.rstrip().split()[1])[1] == '-'):
            data[1002][loc] = 0
        elif (str(line.rstrip().split()[1])[0] == str(
                line.rstrip().split()[1])[1]) and (str(
                    line.rstrip().split()[1])[0] in [
                        'A', 'T', 'G', 'C'
                    ]) and (str(
                        line.rstrip().split()[1])[1] in ['A', 'T', 'G', 'C']):
            if manifest.snps[loc].find(str(line.rstrip().split()[1])[0]) != -1:
                data[1002][loc] = manifest.snps[loc].find(
                    str(line.rstrip().split()[1])[0])
            else:
                print(
                    'WARNING! {} allele possibilities do not match manifest.  {}={} and manifest={}. This snp will not be updated.'
                    .format(line.rstrip().split()[0],
                            line.rstrip().split()[0], originalSnp,
                            manifest.snps[loc]))
                sys.stdout.flush()
                data[1003][loc] = originalSnp
        else:
            pass

        return data

    def validateUpdate(originalGTC, outputName, outDir):
        original_genotype = GenotypeCalls(originalGTC)
        gtc_copy = GenotypeCalls(os.path.join(outDir,
                                              '{}.gtc'.format(outputName)),
                                 check_write_complete=False)
        try:
            assert gtc_copy.get_autocall_date(
            ) == original_genotype.get_autocall_date()
            assert gtc_copy.get_autocall_version(
            ) == original_genotype.get_autocall_version()
            #assert gtc_copy.get_base_calls() == genotype_calls.get_base_calls()
            assert gtc_copy.get_call_rate() == original_genotype.get_call_rate(
            )
            assert gtc_copy.get_cluster_file(
            ) == original_genotype.get_cluster_file()
            assert (gtc_copy.get_control_x_intensities() ==
                    original_genotype.get_control_x_intensities()).all()
            assert (gtc_copy.get_control_y_intensities() ==
                    original_genotype.get_control_y_intensities()).all()
            assert gtc_copy.get_num_no_calls(
            ) == original_genotype.get_num_no_calls()
            assert gtc_copy.get_gender() == original_genotype.get_gender()
            assert (gtc_copy.get_genotype_scores() ==
                    original_genotype.get_genotype_scores()).all()
            #assert gtc_copy.get_genotypes() == genotype_calls.get_genotypes()
            assert gtc_copy.get_percentiles_x(
            ) == original_genotype.get_percentiles_x()
            assert (gtc_copy.get_raw_x_intensities() ==
                    original_genotype.get_raw_x_intensities()).all()

            all_genotypes = gtc_copy.get_genotypes()

            assert len(manifest.names) == len(all_genotypes)
            assert len(manifest.names) == len(gtc_copy.get_logr_ratios())
            assert len(manifest.names) == len(gtc_copy.get_ballele_freqs())
            print(
                os.path.join(outDir, '{}.gtc'.format(outputName)) +
                ' passed validation!')
            sys.stdout.flush()

        except AssertionError:
            print(
                os.path.join(outDir, '{}.gtc'.format(outputName)) +
                ' failed validation -- please re-run this gtc')
            sys.stdout.flush()

    manifest = BeadPoolManifest(bpm)
    manifest.snps[manifest.names.index(
        'rs12248560.1')] = '[T/C]'  # known mistake in bpm
    gtc = ''
    total_gtcs = 0
    data = {}

    with open(snpsToUpdate) as updates:
        for line in updates:
            if line[0] == ">":
                total_gtcs += 1
                if total_gtcs == 1:
                    gtc = line.rstrip().split()[0][1:]
                    outputName = line.rstrip().split()[1]
                    data = getGtcInfo(gtc=os.path.join(gtcDir, gtc))
                    if len(line.rstrip().split()
                           ) == 3:  # means there is metadata to update
                        print("Metadata found. Updating metadata...")
                        sys.stdout.flush()
                        data = updateMetaData(
                            data=data, metaData=line.rstrip().split()[2])
                else:
                    print('Writing updated GTC to new GTC file...')
                    sys.stdout.flush()
                    with open(
                            os.path.join(outDir, '{}.gtc'.format(outputName)),
                            "wb") as output_handle:
                        write_gtc.write_gtc(data, output_handle)

                    validateUpdate(originalGTC=os.path.join(gtcDir, gtc),
                                   outDir=outDir,
                                   outputName=outputName)

                    gtc = line.rstrip().split()[0][1:]
                    outputName = line.rstrip().split()[1]
                    data = getGtcInfo(gtc=os.path.join(gtcDir, gtc))
                    if len(line.rstrip().split()
                           ) == 3:  # means there is metadata to update
                        print("Metadata found. Updating metadata...")
                        sys.stdout.flush()
                        data = updateMetaData(
                            data=data, metaData=line.rstrip().split()[2])

            else:
                data = snpUpdate(data=data, line=line)

    # always the last gtc because out of lines in file at this point
    print('Writing final updated GTC to new GTC file...')
    sys.stdout.flush()
    with open(os.path.join(outDir, '{}.gtc'.format(outputName)),
              "wb") as output_handle:
        write_gtc.write_gtc(data, output_handle)

    validateUpdate(originalGTC=os.path.join(gtcDir, gtc),
                   outDir=outDir,
                   outputName=outputName)

    print("All processing is finished!")
    sys.exit()


def getSampleInfo(bpm, gtcDir, outDir):
    def getGTCinfo(gtc):
        gtcData = GenotypeCalls(gtc)

        data = {}
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_NAME] = gtcData.get_sample_name(
             )  # key:10
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_PLATE] = gtcData.get_sample_plate()
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_WELL] = gtcData.get_sample_well(
             )
            # key:12
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CALL_RATE] = gtcData.get_call_rate(
            )  # key:1006
        data[GenotypeCalls._GenotypeCalls__ID_GC10] = gtcData.get_gc10(
        )  # key:1009
        data[GenotypeCalls.
             _GenotypeCalls__ID_GENDER] = gtcData.get_gender(
             )  # key:1007
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_LOGR_DEV] = gtcData.get_logr_dev(
            )  # key:1008 
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SNP_MANIFEST] = gtcData.get_snp_manifest(
            )  # key:101 
        return data

    nameMatch = open(os.path.join(outDir, 'sentrixMap.txt'), 'w')
    manifest = BeadPoolManifest(bpm)
    input_gtc_list = [
        gtc for gtc in os.listdir(gtcDir) if gtc.endswith(".gtc")
    ]

    header = ['BTID', 'plate', 'well', 'gtcName', 'sampleID', 'callRate', 'gc10', 'sex', 'logrDev']
    nameMatch.write('\t'.join(header) + '\n')
    
    for sampleGtc in input_gtc_list:
        try:
            names = getGTCinfo(gtc=os.path.join(gtcDir, sampleGtc))
            assert manifest.manifest_name == names[101]
            nameMatch.write(names[10] + '\t' + names[11] + '\t' + names[12] +
                            '\t' + sampleGtc + '\t' +
                            '{}-{}-{}'.format(names[11], names[12], names[10]) + 
                            '\t' + str(names[1006]) + '\t' + str(names[1009]) + '\t' +
                            str(names[1007]) + '\t' + str(names[1008]) +
                            '\n')

        except AssertionError:
            print("Error, sample {} in gtc {} does not have matching manifest/bpm file. Sample manifest is listed as {}.  Skipping sample.".format(
             names[10], sampleGtc, names[101]))

    nameMatch.flush()
    nameMatch.close()


def getControlsIntensity(gtcDir, bpm, outDir):
    
    def getGtcInfo(gtc):
        data = {}
        genotype_calls = GenotypeCalls(os.path.join(gtcDir, gtc))
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CLUSTER_FILE] = genotype_calls.get_cluster_file(
            )  # key:100
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SNP_MANIFEST] = genotype_calls.get_snp_manifest(
            )  # key:101
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_NAME] = genotype_calls.get_sample_name(
             )  # key:10
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CONTROLS_X] = genotype_calls.get_control_x_intensities(
            )  # key:500 - 92 values
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CONTROLS_Y] = genotype_calls.get_control_y_intensities(
            )  # key:501 - 92 values
	data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_WELL] = genotype_calls.get_sample_well(
             )  # key:12
	data[
            GenotypeCalls.
            _GenotypeCalls__ID_SAMPLE_PLATE] = genotype_calls.get_sample_plate(
            )  # key:11


        return data

    manifest = BeadPoolManifest(bpm)


    input_gtc_list = [
        gtc for gtc in os.listdir(gtcDir) if gtc.endswith(".gtc")
    ]

    #TODO: detect which probes are available based on CSV input
    #TODO: use the names in CSV and then get a BCP mapping



    intensity_probes_X = [
        "STAINING_DNP_HIGH_1X", "STAINING_DNP_BGND_1X",
        "STAINING_BIOTIN_HIGH_1X", "STAINING_BIOTIN_BGND_1X", "EXTENSION_A_1X",
        "EXTENSION_T_1X", "EXTENSION_C_1X", "EXTENSION_G_1X",
        "TARGET_REMOVAL_1X", "HYBRIDIZATION_HYB_HIGH_1X",
        "HYBRIDIZATION_HYB_MEDIUM_1X", "HYBRIDIZATION_HYB_LOW_1X",
        "STRINGENCY_STRING_PM_1X", "STRINGENCY_STRING_MM_1X",
        "NSB_BGND_RED_1X", "NSB_BGNF_PURPLE_1X", "NSB_BGND_BLUE_1X",
        "NSB_BGND_GREEN_1X", "NON_POLYMORPHIC_NP_A_1X",
        "NON_POLYMORPHIC_NP_T_1X", "NON_POLYMORPHIC_NP_C_1X",
        "NON_POLYMORPHIC_NP_G_1X", "RESTORE_X"
    ]

    intensity_probes_Y = [
        "STAINING_DNP_HIGH_1Y", "STAINING_DNP_BGND_1Y",
        "STAINING_BIOTIN_HIGH_1Y", "STAINING_BIOTIN_BGND_1Y", "EXTENSION_A_1Y",
        "EXTENSION_T_1Y", "EXTENSION_C_1Y", "EXTENSION_G_1Y",
        "TARGET_REMOVAL_1Y", "HYBRIDIZATION_HYB_HIGH_1Y",
        "HYBRIDIZATION_HYB_MEDIUM_1Y", "HYBRIDIZATION_HYB_LOW_1Y",
        "STRINGENCY_STRING_PM_1Y", "STRINGENCY_STRING_MM_1Y",
        "NSB_BGND_RED_1Y", "NSB_BGNF_PURPLE_1Y", "NSB_BGND_BLUE_1Y",
        "NSB_BGND_GREEN_1Y", "NON_POLYMORPHIC_NP_A_1Y",
        "NON_POLYMORPHIC_NP_T_1Y", "NON_POLYMORPHIC_NP_C_1Y",
        "NON_POLYMORPHIC_NP_G_1Y", "RESTORE_Y"
    ]

    intensities_per_sample = {}

    for gtc in input_gtc_list:
        data = getGtcInfo(gtc=gtc)
        try:
            assert data[101] == manifest.manifest_name
            intensities_per_sample['{}-{}-{}'.format(data[11], data[12], data[10])] = {}
            intensityIndex = 0
            for i in range(0, len(data[500])):
                if i%4 == 0:
                    intensities_per_sample['{}-{}-{}'.format(data[11], data[12], data[10])][intensity_probes_X[intensityIndex]] = data[500][i]
                    intensities_per_sample['{}-{}-{}'.format(data[11], data[12], data[10])][intensity_probes_Y[intensityIndex]] = data[501][i]
                    intensityIndex += 1
                else:
                    continue
        except AssertionError:
            print("Sample {}, {} does not have a matching MEGA2 bpm. Skipping sample.".format(data[10], gtc))
            sys.stdout.flush()
            continue
            

    allIntesities = pandas.DataFrame(intensities_per_sample)
    allIntensities_transpose = allIntesities.T
    if os.path.exists(outDir) == False:
        os.mkdir(outDir)
    allIntensities_transpose.to_csv(os.path.join(outDir, 'intesities_per_sample.txt'), index = True, sep = '\t')

def main(bpm, gtcDir, outDir):

    manifest = BeadPoolManifest(bpm)
    input_gtc_list = [
        gtc for gtc in os.listdir(gtcDir) if gtc.endswith(".gtc")
    ]

    for samples in input_gtc_list:
        genotype_calls = GenotypeCalls(os.path.join(gtcDir, samples))

        data = {}
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_AUTOCALL_DATE] = genotype_calls.get_autocall_date(
            )  # key:201
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_AUTOCALL_VERSION] = genotype_calls.get_autocall_version(
            )  # key:300
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_B_ALLELE_FREQS] = genotype_calls.get_ballele_freqs(
            )  # key:1012
        data[GenotypeCalls.
             _GenotypeCalls__ID_BASE_CALLS] = genotype_calls.get_base_calls(
             )  # key:1003
        data[GenotypeCalls.
             _GenotypeCalls__ID_CALL_RATE] = genotype_calls.get_call_rate(
             )  # key:1006
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CLUSTER_FILE] = genotype_calls.get_cluster_file(
            )  # key:100
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CONTROLS_X] = genotype_calls.get_control_x_intensities(
            )  # key:500
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_CONTROLS_Y] = genotype_calls.get_control_y_intensities(
            )  # key:501
        data[GenotypeCalls._GenotypeCalls__ID_GC10] = genotype_calls.get_gc10(
        )  # key:1009
        data[GenotypeCalls._GenotypeCalls__ID_GC50] = (
            genotype_calls.get_gc50(), genotype_calls.get_num_calls(),
            genotype_calls.get_num_no_calls(),
            genotype_calls.get_num_intensity_only())  # key:1011
        data[GenotypeCalls.
             _GenotypeCalls__ID_GENDER] = genotype_calls.get_gender(
             )  # key:1007
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_GENOTYPE_SCORES] = genotype_calls.get_genotype_scores(
            )  # key:1004
        data[GenotypeCalls.
             _GenotypeCalls__ID_GENOTYPES] = genotype_calls.get_genotypes(
             )  # key:1002
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_IMAGING_DATE] = genotype_calls.get_imaging_date(
            )  # key:200
        data[GenotypeCalls.
             _GenotypeCalls__ID_LOGR_DEV] = genotype_calls.get_logr_dev(
             )  # key:1008
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_NORMALIZATION_TRANSFORMS] = genotype_calls.get_normalization_transforms(
            )  # key:400
        data[GenotypeCalls.
             _GenotypeCalls__ID_NUM_SNPS] = genotype_calls.get_num_snps(
             )  # key:1
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_PERCENTILES_X] = genotype_calls.get_percentiles_x(
            )  # key:1014
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_PERCENTILES_Y] = genotype_calls.get_percentiles_y(
            )  # key:1015
        data[GenotypeCalls.
             _GenotypeCalls__ID_PLOIDY] = genotype_calls.get_ploidy()  # key:2
        data[GenotypeCalls.
             _GenotypeCalls__ID_PLOIDY_TYPE] = genotype_calls.get_ploidy_type(
             )  # key:3
        data[GenotypeCalls.
             _GenotypeCalls__ID_RAW_X] = genotype_calls.get_raw_x_intensities(
             )  # key:1000
        data[GenotypeCalls.
             _GenotypeCalls__ID_RAW_Y] = genotype_calls.get_raw_y_intensities(
             )  # key:1001
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_NAME] = genotype_calls.get_sample_name(
             )  # key:10
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SAMPLE_PLATE] = genotype_calls.get_sample_plate(
            )  # key:11
        data[GenotypeCalls.
             _GenotypeCalls__ID_SAMPLE_WELL] = genotype_calls.get_sample_well(
             )  # key:12
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SCANNER_DATA] = genotype_calls.get_scanner_data(
            )  # key:1005
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SLIDE_IDENTIFIER] = genotype_calls.get_slide_identifier(
            )  # key:1016
        data[
            GenotypeCalls.
            _GenotypeCalls__ID_SNP_MANIFEST] = genotype_calls.get_snp_manifest(
            )  # key:101
        data[GenotypeCalls.
             _GenotypeCalls__ID_LOGR_RATIOS] = genotype_calls.get_logr_ratios(
             )  # key:1013

        all_genotypes = genotype_calls.get_genotypes()

        #print(genotype_calls.get_base_calls())

        assert len(manifest.names) == len(all_genotypes)
        assert len(manifest.names) == len(data[1013])
        assert len(manifest.names) == len(data[1012])

        df = pandas.DataFrame({
            'marker name': manifest.names,
            'Chr': manifest.chroms,
            'pos': manifest.map_infos,
            'b allele freq': data[1012],
            'log r ratio': data[1013],
            'base calls': data[1003],
            'genotype calls': data[1002]
        })

        df.to_csv(os.path.join(outDir, samples[:-4] + '_summary_output.csv'),
                  index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Extracts information from gtc files')
    parser.add_argument(
        '--gtcDir',
        type=str,
        default=os.getcwd(),
        help=
        'Full path to location of directory/folder containing gtc files to process (files must end in .gtc)'
    )
    parser.add_argument(
        '--bpm',
        required=True,
        type=str,
        help=
        'Full path to bead pool manifest file (.bpm); must be same one used to generate gtc'
    )
    parser.add_argument(
        '--outDir',
        default=os.getcwd(),
        type=str,
        help='Full path to directory or folder to output results')
    parser.add_argument('--snpUpdates',
                        default=None,
                        type=str,
                        help="Full path to file containing snps to update")

    args = parser.parse_args()
    '''
    Uncomment commands below to activate
    '''

    #main(bpm=args.bpm, gtcDir=args.gtcDir, outDir=args.outDir)
    #manipulate_gtc(bpm=args.bpm, gtcDir=args.gtcDir, snpsToUpdate=args.snpUpdates, outDir=args.outDir)
    #getSampleInfo(bpm=args.bpm, gtcDir=args.gtcDir, outDir=args.outDir)
    getControlsIntensity(gtcDir=args.gtcDir, bpm=args.bpm, outDir=args.outDir)
