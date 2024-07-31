
def extract_clipped_reads(sam, window_info):
    """
    | EXTRACT BREAKPOINT-SUPPORTING CLIPPED READS
    : sam : sam file object of pysam
    : contig : variant contig
    : start : window start position
    : end : window end position
    """
    contig, mate_contig, breakpoint, matepoint, start, end = window_info
    clip_l = set()
    for read in sam.fetch(contig, start, end):
        if ((not read.is_unmapped) and (read.mapping_quality >= 20) and (read.mapping_quality != 255) \
        and (not read.flag & 0x100) \
        and (not read.flag & 0x400)): # 0x100 : not primary-alignment(secondary alignment), 0x400 : read is PCR or optical duplicate
            cigar = read.cigarstring
            align_start = read.reference_start
            align_end = read.reference_end -1
            if 'S' in cigar:
                cigar_match = re.findall(r'(\d+)([A-Z])', cigar)
                clip_findex = cigar.find('S')
                clip_rindex = cigar.rfind('S')
                match_index = cigar.find('M')

                n_match = [int(match[0]) for match in cigar_match if match[1] not in ('S', 'I')]
                len_match = sum(n_match)

                if (clip_findex > match_index) and (clip_rindex > match_index):
                    clip_pos = align_start + len_match -1

                elif (clip_findex < match_index) and (clip_rindex > match_index):
                    if breakpoint < matepoint:
                        clip_pos = align_start
                    else:
                        clip_pos = align_start + len_match -1

                elif (clip_findex < match_index) and (clip_rindex < match_index):
                    clip_pos = align_start

                if clip_pos == breakpoint:
                    #print(contig, clip_pos, breakpoint, read.query_name)
                    clip_l.add(read.query_name)
    return clip_l



def get_clipped_read_count(bam, vcf):
    """
    | GET BREAKPOINT-SUPPORTING CLIPPED READ COUNTS
    : bam : input bam file
    : vcf : input vcf file
    """
    sam_obj = pysam.AlignmentFile(bam, 'rb')
    v_parser = md.VariantParser(vcf)
    v_parser.parse()

    v_l = v_parser.get_variant_list()
    header_l = v_parser.get_header()

    for v in v_l:
        vwindow_l = md.get_variant_window(v, win=150)
        v.set_variant_window(vwindow_l)

        for i, w in enumerate(vwindow_l):
            contig, mate_contig, breakpoint, matepoint, wstart, wend = w
            clip_l = extract_clipped_reads(sam_obj, contig, wstart, wend, breakpoint)
            v.vsplit_l.append(len(clip_l))

    sam_obj.close()
    return v_l, header_l


def filter_by_clipped_read_count(variant_l, child_clip_l, maternal_clip_l, paternal_clip_l):

    pass
    return None

def get_trio_clipped_read_count(child, maternal, paternal, vcf):
    """
    | GET BREAKPOINT-SUPPORTING CLIPPED READ COUNTS
    : bam : input bam file
    : vcf : input vcf file
    """
    child_sam = pysam.AlignmentFile(child, 'rb')
    maternal_sam = pysam.AlignmentFile(maternal, 'rb')
    paternal_sam = pysam.AlignmentFile(paternal, 'rb')

#    sam_obj = pysam.AlignmentFile(bam, 'rb')
    v_parser = md.VariantParser(vcf)
    v_parser.parse()

    v_l = v_parser.get_variant_list()
    header_l = v_parser.get_header()

    for v in v_l:
        vwindow_l = md.get_variant_window(v, win=150)
        v.set_variant_window(vwindow_l)

        for i, w in enumerate(vwindow_l):
            child_clip_l = extract_clipped_reads(child_sam, w)
            maternal_clip_l = extract_clipped_reads(maternal_sam, w)
            paternal_clip_l = extract_clipped_reads(paternal_sam, w)
#            child_clip_l, matenral_clip_l, paternal_clip_l = get_trio_clipped_read_count(child_sam, maternal_sam, paternal_sam, w)
            if len(maternal_clip_l) < 3 and len(paternal_clip_l) < 3:
                v.germline_flag_l.append(False)
#                v.vsplit_l.append(len(child_clip_l))
#                v.msplit_l.append(len(maternal_clip_l))
#                v.psplit_l.append(len(paternal_clip_l))
            else:
                filter_flag = True
                v.germline_flag_l.append(True)
                #break
            v.vsplit_l.append(len(child_clip_l))
            v.msplit_l.append(len(maternal_clip_l))
            v.psplit_l.append(len(paternal_clip_l))
   
     
    child_sam.close()
    maternal_sam.close()
    paternal_sam.close()

    return v_l, header_l



def write_vcf(v_l, header_l, output, total_output):
    """
    | WRITE AN OUTPUT VCF FILE WITH CLIPPED READ COUNT ADDED IN THE INFO FIELD
    : v_l : list containing variant objects
    : header_l : a list containing the header lines of the input vcf file
    : output : output vcf file with clipped-read counts added in the INFO field as CR1, CR2
    """

    filtered_output = output.split('.vcf')
    with open(output, 'w') as outfile, open(total_output, 'w') as total_outfile:
        for h in header_l:
            outfile.write(h)
            total_outfile.write(h)

        for v in v_l:
#            if True in v.germline_flag_l:   continue
            v.add_info('cCR1', v.vsplit_l[0])
            v.add_info('cCR2', v.vsplit_l[1])
            v.add_info('mCR1', v.msplit_l[0])
            v.add_info('mCR2', v.msplit_l[1])
            v.add_info('pCR1', v.psplit_l[0])
            v.add_info('pCR2', v.psplit_l[1])

            line = v.make_vcf_line()
            total_outfile.write(f'{line}\n')
            if not True in v.germline_flag_l:
                outfile.write(f'{line}\n')


def main(args):
    v_l, header_l = get_trio_clipped_read_count(args.bam, args.maternal, args.paternal, args.vcf)
    total_output = args.vcf.split('.vcf')[0] + '.addClipRead.vcf'
    output = total_output.split('.vcf')[0] + '.clipRead_filtered.vcf'
    write_vcf(v_l, header_l, output, total_output)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for processing VCF and calculating performance ")
    parser.add_argument("-i", "--bam", help = " input child bam file ", required = True)
    parser.add_argument("-m", "--maternal", help = " input maternal bam file ", required = True)
    parser.add_argument("-p", "--paternal", help = " input paternal bam file ", required = True)
    parser.add_argument("-v", "--vcf", help = " input vcf file ", required = True)
#    parser.add_argument("-o", "--output", help = " output table ", required = True)
    args = parser.parse_args()

    import re, pysam
    import dnsv_module as md

    main(args)
