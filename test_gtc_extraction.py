from IlluminaBeadArrayFiles import *
import struct
from io import BytesIO
import os
import sys
import argparse
import pandas
import write_gtc


def manipulate_gtc(bpm, snpsToUpdate, outDir):    
    
    def getGtcInfo(gtc): 
        genotype_calls = GenotypeCalls(gtc)
        data[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_DATE] = genotype_calls.get_autocall_date()  # key:201
        data[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_VERSION] = genotype_calls.get_autocall_version()  # key:300
        data[GenotypeCalls._GenotypeCalls__ID_B_ALLELE_FREQS] = genotype_calls.get_ballele_freqs()  # key:1012
        data[GenotypeCalls._GenotypeCalls__ID_BASE_CALLS] = genotype_calls.get_base_calls()  # key:1003
        data[GenotypeCalls._GenotypeCalls__ID_CALL_RATE] = genotype_calls.get_call_rate()  # key:1006
        data[GenotypeCalls._GenotypeCalls__ID_CLUSTER_FILE] = genotype_calls.get_cluster_file()  # key:100
        data[GenotypeCalls._GenotypeCalls__ID_CONTROLS_X] = genotype_calls.get_control_x_intensities()  # key:500
        data[GenotypeCalls._GenotypeCalls__ID_CONTROLS_Y] = genotype_calls.get_control_y_intensities()  # key:501
        data[GenotypeCalls._GenotypeCalls__ID_GC10] = genotype_calls.get_gc10()  # key:1009
        data[GenotypeCalls._GenotypeCalls__ID_GC50] = (genotype_calls.get_gc50(), genotype_calls.get_num_calls(),genotype_calls.get_num_no_calls(),genotype_calls.get_num_intensity_only())  # key:1011
        data[GenotypeCalls._GenotypeCalls__ID_GENDER] = genotype_calls.get_gender()  # key:1007
        data[GenotypeCalls._GenotypeCalls__ID_GENOTYPE_SCORES] = genotype_calls.get_genotype_scores()  # key:1004
        data[GenotypeCalls._GenotypeCalls__ID_GENOTYPES] = genotype_calls.get_genotypes()  # key:1002
        data[GenotypeCalls._GenotypeCalls__ID_IMAGING_DATE] = genotype_calls.get_imaging_date()  # key:200
        data[GenotypeCalls._GenotypeCalls__ID_LOGR_DEV] = genotype_calls.get_logr_dev()  # key:1008
        data[GenotypeCalls._GenotypeCalls__ID_NORMALIZATION_TRANSFORMS] = genotype_calls.get_normalization_transforms()  # key:400
        data[GenotypeCalls._GenotypeCalls__ID_NUM_SNPS] = genotype_calls.get_num_snps()  # key:1
        data[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_X] = genotype_calls.get_percentiles_x()  # key:1014
        data[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_Y] = genotype_calls.get_percentiles_y()  # key:1015
        data[GenotypeCalls._GenotypeCalls__ID_PLOIDY] = genotype_calls.get_ploidy()  # key:2
        data[GenotypeCalls._GenotypeCalls__ID_PLOIDY_TYPE] = genotype_calls.get_ploidy_type()  # key:3
        data[GenotypeCalls._GenotypeCalls__ID_RAW_X] = genotype_calls.get_raw_x_intensities()  # key:1000
        data[GenotypeCalls._GenotypeCalls__ID_RAW_Y] = genotype_calls.get_raw_y_intensities()  # key:1001
        data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_NAME] = genotype_calls.get_sample_name()  # key:10
        data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_PLATE] = genotype_calls.get_sample_plate()  # key:11
        data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_WELL] = genotype_calls.get_sample_well()  # key:12
        data[GenotypeCalls._GenotypeCalls__ID_SCANNER_DATA] = genotype_calls.get_scanner_data()  # key:1005
        data[GenotypeCalls._GenotypeCalls__ID_SLIDE_IDENTIFIER] = genotype_calls.get_slide_identifier()  # key:1016
        data[GenotypeCalls._GenotypeCalls__ID_SNP_MANIFEST] = genotype_calls.get_snp_manifest()  # key:101
        data[GenotypeCalls._GenotypeCalls__ID_LOGR_RATIOS] = genotype_calls.get_logr_ratios()  # key:1013

        return data

    def snpUpdate(data, line):
        loc = manifest.names.index(line.rstrip().split()[0])
        data[1003][loc] = str(line.rstrip().split()[1])
        if ((str(line.rstrip().split()[1])[0] != str(line.rstrip().split()[1])[1]) and (str(line.rstrip().split()[1])[0] != '-')):
            data[1002][loc] = 2
        elif (str(line.rstrip().split()[1])[0] == '-') and (str(line.rstrip().split()[1])[1] == '-'):
            data[1002][loc] = 0
        elif (str(line.rstrip().split()[1])[0] == str(line.rstrip().split()[1])[1]) and (str(line.rstrip().split()[1])[0] in ['A', 'T', 'G', 'C']) and (str(line.rstrip().split()[1])[1] in ['A', 'T', 'G', 'C']):
            data[1002][loc] = manifest.snps[loc].find(str(line.rstrip().split()[1])[0])
        else:
            pass

        return data

    def validateUpdate(originalGTC, outputName, outDir):
        original_genotype = GenotypeCalls(originalGTC) 
        gtc_copy = GenotypeCalls(os.path.join(outDir, '{}.gtc'.format(outputName)), check_write_complete=False)
        try:
            assert gtc_copy.get_autocall_date() == original_genotype.get_autocall_date()
            assert gtc_copy.get_autocall_version() == original_genotype.get_autocall_version()
            #assert gtc_copy.get_base_calls() == genotype_calls.get_base_calls()
            assert gtc_copy.get_call_rate() == original_genotype.get_call_rate()
            assert gtc_copy.get_cluster_file() == original_genotype.get_cluster_file()
            assert (gtc_copy.get_control_x_intensities() == original_genotype.get_control_x_intensities()).all()
            assert (gtc_copy.get_control_y_intensities() == original_genotype.get_control_y_intensities()).all()
            assert gtc_copy.get_num_no_calls() == original_genotype.get_num_no_calls()
            assert gtc_copy.get_gender() == original_genotype.get_gender()
            assert (gtc_copy.get_genotype_scores() == original_genotype.get_genotype_scores()).all()
            #assert gtc_copy.get_genotypes() == genotype_calls.get_genotypes()
            assert gtc_copy.get_percentiles_x() == original_genotype.get_percentiles_x()
            assert (gtc_copy.get_raw_x_intensities() == original_genotype.get_raw_x_intensities()).all()

            all_genotypes = gtc_copy.get_genotypes()
            
            assert len(manifest.names) == len(all_genotypes)
            assert len(manifest.names) == len(gtc_copy.get_logr_ratios())
            assert len(manifest.names) == len(gtc_copy.get_ballele_freqs())
            print(os.path.join(outDir, '{}.gtc'.format(outputName)) + ' passed validation!')
            sys.stdout.flush()

        except AssertionError:
            print(os.path.join(outDir, '{}.gtc'.format(outputName)) + ' failed validation -- please re-run this gtc')
            sys.stdout.flush()

    
    manifest = BeadPoolManifest(bpm)
    manifest.snps[manifest.names.index('rs12248560.1')] = '[T/C]' # known mistake in bpm
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
                    data = getGtcInfo(gtc=gtc)
                else:
                    print('Writing updated GTC to new GTC file...')
                    sys.stdout.flush()
                    with open(os.path.join(outDir, '{}.gtc'.format(outputName)), "wb") as output_handle:
                        write_gtc.write_gtc(data, output_handle)
                    
                    validateUpdate(originalGTC=gtc, outDir=outDir, outputName=outputName)

                    gtc = line.rstrip().split()[0][1:]
                    outputName = line.rstrip().split()[1]
                    data = getGtcInfo(gtc=gtc)
            else:
                data = snpUpdate(data=data, line=line)

    # always the last gtc because out of lines in file at this point    
    print('Writing final updated GTC to new GTC file...')
    sys.stdout.flush()
    with open(os.path.join(outDir, '{}.gtc'.format(outputName)), "wb") as output_handle:
        write_gtc.write_gtc(data, output_handle)
    
    validateUpdate(originalGTC=gtc, outDir=outDir, outputName=outputName)

    print("All processing is finished!")
    sys.exit()

def main(bpm, gtcDir, outDir):

    manifest = BeadPoolManifest(bpm)
    input_gtc_list = [gtc for gtc in os.listdir(gtcDir) if gtc.endswith(".gtc")]
   
    for samples in input_gtc_list:
        genotype_calls = GenotypeCalls(os.path.join(gtcDir, samples))

        data = {}
        data[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_DATE] = genotype_calls.get_autocall_date() # key:201
        data[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_VERSION] = genotype_calls.get_autocall_version()  # key:300
        data[GenotypeCalls._GenotypeCalls__ID_B_ALLELE_FREQS] = genotype_calls.get_ballele_freqs()  # key:1012
        data[GenotypeCalls._GenotypeCalls__ID_BASE_CALLS] = genotype_calls.get_base_calls()  # key:1003
        data[GenotypeCalls._GenotypeCalls__ID_CALL_RATE] = genotype_calls.get_call_rate()  # key:1006
        data[GenotypeCalls._GenotypeCalls__ID_CLUSTER_FILE] = genotype_calls.get_cluster_file()  # key:100
        data[GenotypeCalls._GenotypeCalls__ID_CONTROLS_X] = genotype_calls.get_control_x_intensities()  # key:500
        data[GenotypeCalls._GenotypeCalls__ID_CONTROLS_Y] = genotype_calls.get_control_y_intensities()  # key:501
        data[GenotypeCalls._GenotypeCalls__ID_GC10] = genotype_calls.get_gc10()  # key:1009
        data[GenotypeCalls._GenotypeCalls__ID_GC50] = (genotype_calls.get_gc50(), genotype_calls.get_num_calls(),genotype_calls.get_num_no_calls(),genotype_calls.get_num_intensity_only())  # key:1011
        data[GenotypeCalls._GenotypeCalls__ID_GENDER] = genotype_calls.get_gender()  # key:1007
        data[GenotypeCalls._GenotypeCalls__ID_GENOTYPE_SCORES] = genotype_calls.get_genotype_scores()  # key:1004
        data[GenotypeCalls._GenotypeCalls__ID_GENOTYPES] = genotype_calls.get_genotypes()  # key:1002
        data[GenotypeCalls._GenotypeCalls__ID_IMAGING_DATE] = genotype_calls.get_imaging_date()  # key:200
        data[GenotypeCalls._GenotypeCalls__ID_LOGR_DEV] = genotype_calls.get_logr_dev()  # key:1008
        data[GenotypeCalls._GenotypeCalls__ID_NORMALIZATION_TRANSFORMS] = genotype_calls.get_normalization_transforms()  # key:400
        data[GenotypeCalls._GenotypeCalls__ID_NUM_SNPS] = genotype_calls.get_num_snps()  # key:1
        data[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_X] = genotype_calls.get_percentiles_x()  # key:1014
        data[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_Y] = genotype_calls.get_percentiles_y()  # key:1015
        data[GenotypeCalls._GenotypeCalls__ID_PLOIDY] = genotype_calls.get_ploidy()  # key:2
        data[GenotypeCalls._GenotypeCalls__ID_PLOIDY_TYPE] = genotype_calls.get_ploidy_type()  # key:3
        data[GenotypeCalls._GenotypeCalls__ID_RAW_X] = genotype_calls.get_raw_x_intensities()  # key:1000
        data[GenotypeCalls._GenotypeCalls__ID_RAW_Y] = genotype_calls.get_raw_y_intensities()  # key:1001
        data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_NAME] = genotype_calls.get_sample_name()  # key:10
        data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_PLATE] = genotype_calls.get_sample_plate()  # key:11
        data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_WELL] = genotype_calls.get_sample_well()  # key:12
        data[GenotypeCalls._GenotypeCalls__ID_SCANNER_DATA] = genotype_calls.get_scanner_data()  # key:1005
        data[GenotypeCalls._GenotypeCalls__ID_SLIDE_IDENTIFIER] = genotype_calls.get_slide_identifier()  # key:1016
        data[GenotypeCalls._GenotypeCalls__ID_SNP_MANIFEST] = genotype_calls.get_snp_manifest()  # key:101
        data[GenotypeCalls._GenotypeCalls__ID_LOGR_RATIOS] = genotype_calls.get_logr_ratios()  # key:1013

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

    parser = argparse.ArgumentParser(description='Extracts information from gtc files')
    parser.add_argument('--gtcDir', type=str, help='Full path to location of directory/folder containing gtc files to process (files must end in .gtc)')
    parser.add_argument('--bpm', required=True, type=str, help='Full path to bead pool manifest file (.bpm); must be same one used to generate gtc')
    parser.add_argument('--outDir', default=os.getcwd(), type=str, help='Full path to directory or folder to output results')
    parser.add_argument('--snpUpdates', default=None, type=str, help="Full path to file containing snps to update")

    args = parser.parse_args()

    #main(bpm=args.bpm, gtcDir=args.gtcDir, outDir=args.outDir)
    manipulate_gtc(bpm=args.bpm, snpsToUpdate=args.snpUpdates, outDir=args.outDir)