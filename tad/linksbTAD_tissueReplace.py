import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('target_path')

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

    tissues = {}
    with open('tissuetable_10092018.txt', 'r') as tissue_file:
        for line in tissue_file:
            line_split = line.strip().split('\t')
            tissueID = line_split[0]
            tissue = line_split[1]
            tissue = tissue.replace(' ', '_')
            tissues[tissue] = tissueID

    tissuesReplace(os.path.join(target_path, 'linksbTAD'), os.path.join(target_path, 'linksbTADtissues'))