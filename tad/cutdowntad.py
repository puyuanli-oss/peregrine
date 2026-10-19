import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('target_path')

def cutdowntad(eqtl_file, heirarchical_tad_file, tad_file, output_file):
    '''
        only want to include the TAD links that support existing links from ChIA-PET or eQTL data
    '''
    hash = {}
    count1 = 0
    count2 = 0

    def line_parse(line, hash):
        line_split = line.strip().split('\t')
        enhancer = line_split[0]
        panthid = line_split[1]
        link = enhancer + '_' + panthid
        hash[link] = 1

    with open(eqtl_file, 'r') as eqtl:
        for line in eqtl:
            line_parse(line, hash)

    with open(heirarchical_tad_file, 'r') as heirarchical:
        for line in heirarchical:
            line_parse(line, hash)

    with open(output_file, 'w') as out, open(tad_file, 'r') as tad:
        for line in tad:
            enhancer, panthid, tissue, assay = line.strip().split('\t')
            link = enhancer + '_' + panthid
            rest = f"{enhancer}\t{panthid}\t{tissue}\t{assay}"
            if link in hash:
                count1 += 1
                out.write(f"{rest}\n")
            else:
                count2 += 1

    print(f"Number of TAD links found in TADintxn, eQTL, or ChIA: {count1}")

if __name__ == "__main__":

    args = parser.parse_args()
    target_path = args.target_path
    os.makedirs(target_path, exist_ok=True)

    target_path_base = os.path.dirname(os.path.abspath(target_path))
    cutdowntad(target_path_base + '/eQTL/linksDBnumeqtl', target_path_base + '/heirarchicalTAD/PSYCHIClinksDB', os.path.join(target_path, 'linksTADtissues'), os.path.join(target_path, 'selectTAD'))