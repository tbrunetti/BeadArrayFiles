from IlluminaBeadArrayFiles import GenotypeCalls, BeadArrayUtility
import struct
from io import BytesIO

def handle_int(value):
    return struct.pack("<i", value)

def handle_short(value):
    return struct.pack("<H", value)

def handle_char(value):
    return struct.pack("c", value)

def handle_byte(value):
    #print(value)
    return struct.pack("B", value)

def handle_float(value):
    return struct.pack("<f", value)

def handle_gc50(value):
    return struct.pack("<fiii", value[0], value[1], value[2], value[3])

def handle_percentiles(value):
    return struct.pack("<HHH", value[0], value[1], value[2])

def handle_string(value):
    assert len(value) <= 127
    return struct.pack("B", len(value)) + value

def handle_basecalls(value):
    return value

def handle_scanner_data(value):
    return handle_string(value.name) + handle_int(value.pmt_green) + handle_int(value.pmt_red) + handle_string(value.version) + handle_string(value.user)

def handle_normalization_transform(value):
    return struct.pack("<iffffff", value.version, value.offset_x, value.offset_y, value.scale_x, value.scale_y, value.shear, value.theta)


toc2handler = {}
toc2handler[GenotypeCalls._GenotypeCalls__ID_NUM_SNPS] = handle_int
toc2handler[GenotypeCalls._GenotypeCalls__ID_PLOIDY] = handle_int
toc2handler[GenotypeCalls._GenotypeCalls__ID_PLOIDY_TYPE] = handle_int
toc2handler[GenotypeCalls._GenotypeCalls__ID_SAMPLE_NAME] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_SAMPLE_PLATE] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_SAMPLE_WELL] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_CLUSTER_FILE] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_SNP_MANIFEST] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_IMAGING_DATE] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_DATE] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_VERSION] = handle_string
toc2handler[GenotypeCalls._GenotypeCalls__ID_NORMALIZATION_TRANSFORMS] = handle_normalization_transform
toc2handler[GenotypeCalls._GenotypeCalls__ID_CONTROLS_X] = handle_short
toc2handler[GenotypeCalls._GenotypeCalls__ID_CONTROLS_Y] = handle_short
toc2handler[GenotypeCalls._GenotypeCalls__ID_RAW_X] = handle_short
toc2handler[GenotypeCalls._GenotypeCalls__ID_RAW_Y] = handle_short
toc2handler[GenotypeCalls._GenotypeCalls__ID_GENOTYPES] = handle_byte
toc2handler[GenotypeCalls._GenotypeCalls__ID_BASE_CALLS] = handle_basecalls
toc2handler[GenotypeCalls._GenotypeCalls__ID_GENOTYPE_SCORES] = handle_float
toc2handler[GenotypeCalls._GenotypeCalls__ID_SCANNER_DATA] = handle_scanner_data
toc2handler[GenotypeCalls._GenotypeCalls__ID_CALL_RATE] = handle_float
toc2handler[GenotypeCalls._GenotypeCalls__ID_GENDER] = handle_char
toc2handler[GenotypeCalls._GenotypeCalls__ID_LOGR_DEV] = handle_float
toc2handler[GenotypeCalls._GenotypeCalls__ID_GC10] = handle_float
toc2handler[GenotypeCalls._GenotypeCalls__ID_GC50] = handle_gc50
toc2handler[GenotypeCalls._GenotypeCalls__ID_B_ALLELE_FREQS] = handle_float
toc2handler[GenotypeCalls._GenotypeCalls__ID_LOGR_RATIOS] = handle_float
toc2handler[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_X] = handle_percentiles
toc2handler[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_Y] = handle_percentiles
toc2handler[GenotypeCalls._GenotypeCalls__ID_SLIDE_IDENTIFIER] = handle_string

list_types = []
list_types.append(GenotypeCalls._GenotypeCalls__ID_NORMALIZATION_TRANSFORMS)
list_types.append(GenotypeCalls._GenotypeCalls__ID_CONTROLS_X)
list_types.append(GenotypeCalls._GenotypeCalls__ID_CONTROLS_Y)
list_types.append(GenotypeCalls._GenotypeCalls__ID_RAW_X)
list_types.append(GenotypeCalls._GenotypeCalls__ID_RAW_Y)
list_types.append(GenotypeCalls._GenotypeCalls__ID_GENOTYPES)
list_types.append(GenotypeCalls._GenotypeCalls__ID_BASE_CALLS)
list_types.append(GenotypeCalls._GenotypeCalls__ID_GENOTYPE_SCORES)
list_types.append(GenotypeCalls._GenotypeCalls__ID_B_ALLELE_FREQS)
list_types.append(GenotypeCalls._GenotypeCalls__ID_LOGR_RATIOS)


def write_gtc(data, handle):
    handle.write(b'g')
    handle.write(b't')
    handle.write(b'c')
    handle.write(handle_byte(5))

    num_entries = len(data)
    handle.write(handle_int(num_entries))
    offset = 8 + num_entries * 6;


    buffer = BytesIO()
    for toc_id in data:
        # write the toc ID
        handle.write(handle_short(toc_id))

        # write the data into the buffer
        if toc_id in list_types:
            handle.write(handle_int(offset + buffer.tell()))
            buffer.write(handle_int(len(data[toc_id])))
            for element in data[toc_id]:
                buffer.write(toc2handler[toc_id](element))
        else:
            if toc2handler[toc_id] == handle_int:
                handle.write(handle_int(data[toc_id]))
            else:
                handle.write(handle_int(offset + buffer.tell()))
                buffer.write(toc2handler[toc_id](data[toc_id]))
    buffer.seek(0)
    handle.write(buffer.read1(-1))

'''
def main(input_gtc):
    genotype_calls = GenotypeCalls(input_gtc)

    data = {}
    data[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_DATE] = genotype_calls.get_autocall_date()
    data[GenotypeCalls._GenotypeCalls__ID_AUTOCALL_VERSION] = genotype_calls.get_autocall_version()
    data[GenotypeCalls._GenotypeCalls__ID_B_ALLELE_FREQS] = genotype_calls.get_ballele_freqs()
    data[GenotypeCalls._GenotypeCalls__ID_BASE_CALLS] = genotype_calls.get_base_calls()
    data[GenotypeCalls._GenotypeCalls__ID_CALL_RATE] = genotype_calls.get_call_rate()
    data[GenotypeCalls._GenotypeCalls__ID_CLUSTER_FILE] = genotype_calls.get_cluster_file()
    data[GenotypeCalls._GenotypeCalls__ID_CONTROLS_X] =  genotype_calls.get_control_x_intensities()
    data[GenotypeCalls._GenotypeCalls__ID_CONTROLS_Y] =  genotype_calls.get_control_y_intensities()
    data[GenotypeCalls._GenotypeCalls__ID_GC10] = genotype_calls.get_gc10()
    data[GenotypeCalls._GenotypeCalls__ID_GC50] = (genotype_calls.get_gc50(), genotype_calls.get_num_calls(), genotype_calls.get_num_no_calls(), genotype_calls.get_num_intensity_only())
    data[GenotypeCalls._GenotypeCalls__ID_GENDER] = genotype_calls.get_gender()
    data[GenotypeCalls._GenotypeCalls__ID_GENOTYPE_SCORES] = genotype_calls.get_genotype_scores()
    data[GenotypeCalls._GenotypeCalls__ID_GENOTYPES] = genotype_calls.get_genotypes()
    data[GenotypeCalls._GenotypeCalls__ID_IMAGING_DATE] = genotype_calls.get_imaging_date()
    data[GenotypeCalls._GenotypeCalls__ID_LOGR_DEV] = genotype_calls.get_logr_dev()
    data[GenotypeCalls._GenotypeCalls__ID_NORMALIZATION_TRANSFORMS] = genotype_calls.get_normalization_transforms()
    data[GenotypeCalls._GenotypeCalls__ID_NUM_SNPS] = genotype_calls.get_num_snps()
    data[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_X] = genotype_calls.get_percentiles_x()
    data[GenotypeCalls._GenotypeCalls__ID_PERCENTILES_Y] = genotype_calls.get_percentiles_y()
    data[GenotypeCalls._GenotypeCalls__ID_PLOIDY] = genotype_calls.get_ploidy()
    data[GenotypeCalls._GenotypeCalls__ID_PLOIDY_TYPE] = genotype_calls.get_ploidy_type()
    data[GenotypeCalls._GenotypeCalls__ID_RAW_X] = genotype_calls.get_raw_x_intensities()
    data[GenotypeCalls._GenotypeCalls__ID_RAW_Y] = genotype_calls.get_raw_y_intensities()
    data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_NAME] = genotype_calls.get_sample_name()
    data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_PLATE] = genotype_calls.get_sample_plate()
    data[GenotypeCalls._GenotypeCalls__ID_SAMPLE_WELL] = genotype_calls.get_sample_well()
    data[GenotypeCalls._GenotypeCalls__ID_SCANNER_DATA] = genotype_calls.get_scanner_data()
    data[GenotypeCalls._GenotypeCalls__ID_SLIDE_IDENTIFIER] = genotype_calls.get_slide_identifier()
    data[GenotypeCalls._GenotypeCalls__ID_SNP_MANIFEST] = genotype_calls.get_snp_manifest()
    
    with open("output.gtc", "wb") as output_handle:
        write_gtc(data, output_handle)

    gtc_copy = GenotypeCalls("output.gtc", check_write_complete=False)
    assert gtc_copy.get_autocall_date() == genotype_calls.get_autocall_date()
    assert gtc_copy.get_autocall_version() == genotype_calls.get_autocall_version()
    assert gtc_copy.get_base_calls() == genotype_calls.get_base_calls()
    assert gtc_copy.get_call_rate() == genotype_calls.get_call_rate()
    assert gtc_copy.get_cluster_file() == genotype_calls.get_cluster_file()
    assert (gtc_copy.get_control_x_intensities() == genotype_calls.get_control_x_intensities()).all()
    assert (gtc_copy.get_control_y_intensities() == genotype_calls.get_control_y_intensities()).all()
    assert gtc_copy.get_num_no_calls() == genotype_calls.get_num_no_calls()
    assert gtc_copy.get_gender() == genotype_calls.get_gender()
    assert (gtc_copy.get_genotype_scores() == genotype_calls.get_genotype_scores()).all()
    assert gtc_copy.get_genotypes() == genotype_calls.get_genotypes()
    assert gtc_copy.get_percentiles_x() == genotype_calls.get_percentiles_x()
    assert (gtc_copy.get_raw_x_intensities() == genotype_calls.get_raw_x_intensities()).all()


if __name__ == "__main__":
    import sys
    main(sys.argv[1])
'''







