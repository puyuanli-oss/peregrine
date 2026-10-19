import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('target_path')

from collections import defaultdict

def TADlinks(gene_file, enhancer_tad_file, output_file):
    '''
        Trim down to a file showing just the gene>enhancer>cell>assay
    '''

    tadgenes = defaultdict(dict)
    with open(gene_file, 'r') as tss_file:
        for line in tss_file:
            chr, start, end, geneID, one, two, three, tad, tissue = line.strip().split('\t')
            cell = tissue.split('|')[1]
            tad = tad + '\t' + cell
            tadgenes[tad][geneID] = 1

    tadenh = defaultdict(dict)
    with open(enhancer_tad_file, 'r') as enh_file:
        for line in enh_file:
            chr, start, end, enhID, one, two, three, tad, tissue = line.strip().split('\t')
            cell = tissue.split('|')[1]
            tad = tad + '\t' + cell
            tadenh[enhID][tad] = 1

    with open(output_file, 'w') as out:
        for enh, tads in tadenh.items():
            for tad in tads:
                if tad in tadgenes:
                    for gene in tadgenes[tad]:
                        panthID = panther_mapping.get(gene, '')
                        cell = tad.split('\t')[1]
                        if enh and panthID and cell:
                            out.write(f"{enh}\t{panthID}\t{cell}\t3\n")

def TADlinks_write_to_file(tadenh, tadgenes, output_file, start_idx, end_idx):
    with open(output_file, 'a') as out:
        for enh in tadenh.keys():
            for tad in list(tadenh[enh].keys())[start_idx:end_idx]:
                if tad in tadgenes:
                    for gene in tadgenes[tad].keys():
                        panthID = panther_mapping.get(gene, '')
                        cell = tad.split('\t')[1]
                        if enh and panthID and cell:
                            out.write(f"{enh}\t{panthID}\t{cell}\t3\n")

def tissuesReplace(input_file, output_file):
    '''
        Replace the tissues and cell types with their codes
    '''
    with open(output_file, 'w') as out:
        with open(input_file, 'r') as tss_file:
            for line in tss_file:
                line_split = line.strip().split('\t')
                if len(line_split) == 4:
                    enhID, panthID, tiss, assay = line_split
                    cell = tissues.get(tiss, '')
                    if not cell.strip():
                        raise ValueError(f"{input_file}: missing tissue code for {tiss!r}")
                    out.write(f"{enhID}\t{panthID}\t{cell}\t{assay}\n")
                else:
                    raise ValueError(f"{input_file}: expected 4 tab-separated fields")
                
if __name__ == "__main__":

    args = parser.parse_args()
    target_path = args.target_path
    os.makedirs(target_path, exist_ok=True)

    panther_mapping = {}
    with open('pantherGeneList.txt', 'r') as pantherGene_file:
        for line in pantherGene_file:
            line_split = line.strip().split('\t')
            panth = line_split[0]
            ensg = line_split[1]
            panther_mapping[ensg] = panth

    tissues = {}
    with open('tissuetable_10092018.txt', 'r') as tissue_file:
        for line in tissue_file:
            line_split = line.strip().split('\t')
            tissueID = line_split[0]
            tissue = line_split[1]
            tissue = tissue.replace(' ', '_')
            tissues[tissue] = tissueID

    TADlinks(os.path.join(target_path, 'TSStTAD'), os.path.join(target_path, 'enhancerstTAD'), os.path.join(target_path, 'linkstTAD'))
    tissuesReplace(os.path.join(target_path, 'linkstTAD'), os.path.join(target_path, 'linkstTADtissues'))