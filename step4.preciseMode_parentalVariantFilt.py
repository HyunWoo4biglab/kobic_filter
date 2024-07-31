# A SCRIPT TO FILTER PARENTAL GERMLINE VARIANTS CALLED BY ETCHING - ETCHING-TRIO PRECISE MODE
# 04-11-2024
# WRITTEN BY HYUN WOO KIM

def get_variants(vcf):

    v_parser = md.VariantParser(vcf)
    v_parser.parse()

    v_l = v_parser.get_variant_list()

    return v_l


def filter_variants(child, paternal, maternal, output, offset=10):
    """
    | A FUNCTION TO COMPARE CHILD VARIANTS WITH ITS PARENTAL VARIATNS CALLED BY ETCHING
    : child : child derived vcf
    : paternal : paternal derived vcf
    : maternal : maternal derived vcf
    : output : output vcf file containing only child specific variants
    : offset : number of offset basepairs to consier as match
    """

    child_v_l = get_variants(child)
    paternal_v_l = get_variants(paternal)
    maternal_v_l = get_variants(maternal)

    parental_v_l = maternal_v_l + paternal_v_l

    with open(output, 'w') as outfile:
        #---- WRITE HEADER
        infile = open(child, 'r')
        for l in infile:
            if l.startswith('#'):
                outfile.write(l)
        infile.close()

        #---- COMPARE WITH PARENTAL VARIANTS
        for v in child_v_l:
            match_bool = False
            for p in parental_v_l:
                if md.variant_match(v, p, offset):
                    match_bool = True
                    break
            if match_bool:  continue
            #else:
            line = v.make_vcf_line()
            outfile.write(f'{line}\n')
        



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for comparing an input VCF with truth benchmark VCF giving a10bp window of start & end coordinate ")
    parser.add_argument("-c", "--child", help = " intput child vcf ", required = True)
    parser.add_argument("-p", "--paternal", help = " paternal vcf ", required = True)
    parser.add_argument("-m", "--maternal", help = " maternal vcf ", required = True)
    parser.add_argument("-o", "--output", help = " output vcf file ", required = True)
    args = parser.parse_args()

    import re
    import dnsv_module as md
 
    filter_variants(args.child, args.paternal, args.maternal, args.output)
