import pandas as pd
from pybedtools import BedTool
import os
import shutil
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('target_path')


def idmapping_processing(id_output_file):
    mapping = {}
    with open('heirarchicalTAD/UP000005640_9606.idmapping', 'r') as file:
        for line in file:
            line_split = line.strip().split('\t')
            if line_split[0] not in mapping:
                mapping[line_split[0]] = {}
            if line_split[1] == 'Gene_Name':
                mapping[line_split[0]]['Gene_Name'] = line_split[2]
            elif line_split[1] == 'HGNC':
                mapping[line_split[0]]['HGNC'] = line_split[2]
            elif line_split[1]== 'Ensembl':
                mapping[line_split[0]]['Ensembl'] = line_split[2].split('.')[0]

    df = pd.DataFrame(columns=['HGNC ID', 'symbol', 'UniProt', 'Ensembl gene ID'])
    count = 0
    for key, val in mapping.items():
        if ('Ensembl' in val) and ('HGNC' in val) and ('Gene_Name') in val:
            row = [val['HGNC'], val['Gene_Name'], key, val['Ensembl']]
            df.loc[count] = row
            count += 1
    df.to_csv(id_output_file, sep='\t', index=False)

def nonegvalues(input_file, output_file, out_file):
    with open(input_file, 'r') as bed_file, open(out_file, 'w') as tossed_out, open(output_file, 'w') as out:
        for line in bed_file:
            line_split = line.strip().split('\t')
            start = line_split[1]
            end = line_split[2]
            if int(start) < 0 or int(end) < 0:
                tossed_out.write(line)
            else:
                out.write(line)


def bedtoolIntersect(cred_file, eqtl_file, output_file):
    cred_file_bedtoolFile = BedTool(cred_file)
    eqtl_file_bedtoolFile = BedTool(eqtl_file)

    intersection = cred_file_bedtoolFile.intersect(eqtl_file_bedtoolFile, wa=True, wb=True, f=0.9)
    with open(output_file, 'w') as out:
        for elt in intersection:
            out.write(str(elt))


def HGNC2PANTH(input_file, output_file, hgnc_file, unmatched_file, tissue_name):
    HGNC_mapping = open(hgnc_file, 'r')
    lines = HGNC_mapping.readlines()
    HGNC_mapping.close()

    hash1 = {}
    hash2 = {}
    hash3 = {}
    for line in lines[1:]:
        hgnc, symbol, uniprot, ensg,  = line.strip().split('\t')
        hgncid = hgnc.split(':')[1]
        hgncid = 'HGNC=' + hgncid
        hash1[symbol] = hgncid
        hash2[symbol] = uniprot
        hash3[symbol] = ensg


    with open(input_file, 'r') as input, open(output_file, 'w') as out, open(unmatched_file, 'w') as unmatched:
        for line in input:
            line = line.rstrip('\n')
            echr, estart, eend, enhID, pchr, pstart, pend, mix = line.strip().split('\t')
            mix_split = mix.split(':')
            gene = mix_split[0]
            pvalue = mix_split[2]

            if gene in hash1:
                for candidate in (hash1[gene], hash2[gene], hash3[gene]):
                    if candidate in panth:
                        panthid = panth[candidate]
                        out.write(f'{echr}\t{estart}\t{eend}\t{enhID}\t{pchr}\t{pstart}\t{pend}\t{panthid}\t{pvalue}\t{tissue_name}\n')
                        break
                else:
                    unmatched.write(f'1\t{line}\n')
            else:
                unmatched.write(f'2\t{line}\n')


def reformat(input_file, output_file):
    hash = {}
    with open(input_file, 'r') as input:
        for line in input:
            echr, estart, eend, enhID, pchr, pstart, pend, gene, pval, tissue = line.strip().split('\t')
            break_ = f"{echr}\t{estart}\t{eend}\t{enhID}\t{gene}\t{tissue}\tTADinteractions\t{pval}"
            hash[break_] = 1

    with open(output_file, 'w') as out:
        for key in sorted(hash.keys()):
            out.write(key + '\n')


def tissuesReplace(input_file, output_file):
    '''
        Replace the tissues and cell types with their codes
    '''
    with open(input_file, 'r') as tss_file, open(output_file, 'w') as out:
        for line in tss_file:
            line_split = line.strip().split('\t')
            enhID = line_split[3]
            gene = line_split[4]
            tissue = line_split[5]
            tissue_code = tissues.get(tissue, '')
            pval = line_split[-1]
            assay = 4
            out.write(f"{enhID}\t{gene}\t{tissue_code}\t{pval}\t{assay}\n")


def concatenate(input_files, output_file):
    '''
        makes a combined file for all tissues
    '''
    with open(output_file, 'w') as outfile:
        for file_name in input_files:
            infile = open(file_name, 'r')
            outfile.write(infile.read())
            infile.close()


if __name__ == "__main__":

    args = parser.parse_args()
    target_path = args.target_path
    os.makedirs(target_path, exist_ok=True)

    panth = {}
    with open('pantherGeneList.txt', 'r') as panth_mapping:
        for line in panth_mapping:
            line_split = line.strip().split()
            longID = line_split[0]
            ensg = line_split[1]
            longID_split = longID.split('|')
            HGNC = longID_split[1]
            uniprotkb = longID_split[2].split('=')[1]
            panth[HGNC] = longID
            panth[uniprotkb] = longID
            panth[ensg] = longID

    tissues = {}
    with open('tissuetable_10092018.txt', 'r') as tissue_file:
        for line in tissue_file:
            line_split = line.strip().split('\t')
            tissueID = line_split[0]
            tissue = line_split[1]
            tissue = tissue.replace(' ', '_')
            tissues[tissue] = tissueID

    id_output_file = os.path.join(target_path, 'resultsHGNC.txt')
    idmapping_processing(id_output_file)

    bed_files = glob.glob(os.path.join('heirarchicalTAD/data', '*'))
    if not bed_files:
        raise FileNotFoundError(
            'No hierarchical TAD input files found in heirarchicalTAD/data'
        )
    db_files = []

    for bed_file in bed_files:
        bed_file_split = bed_file.split('.')

        folder = os.path.basename(bed_file_split[0])
        if '_' in bed_file_split[0]:
            folder = folder.split('_')[1]
            
        folder = os.path.join(target_path, folder)

        tissue = os.path.basename(folder)
        file = os.path.basename(bed_file)
        if tissue in tissues:
            os.makedirs(folder, exist_ok=True)
            shutil.copy(bed_file, folder)

            nonegvalues(os.path.join(folder, file), os.path.join(folder, 'nn_' + file), os.path.join(folder, 'out'))
            file_name = os.path.splitext(file)[0]
            bedtoolIntersect('CREbedDBenhancers_10092018', os.path.join(folder, 'nn_' + file), os.path.join(folder, 'nn_' + file_name + '_intersect'))
            HGNC2PANTH(os.path.join(folder, 'nn_' + file_name + '_intersect'), os.path.join(folder, file_name + '_out'), id_output_file, os.path.join(folder, 'unmatched'), tissue)
            reformat(os.path.join(folder, file_name + '_out'), os.path.join(folder, 'intTADlinks_' + file_name))
            tissuesReplace(os.path.join(folder, 'intTADlinks_' + file_name), os.path.join(folder,  file_name + '_DB'))
            db_files.append(os.path.join(folder,  file_name + '_DB'))

    if not db_files:
        raise ValueError(
            'No hierarchical TAD input files matched the tissue table'
        )
    concatenate(db_files, os.path.join(target_path, 'PSYCHIClinksDB'))
