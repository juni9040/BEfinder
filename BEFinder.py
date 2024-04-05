import numpy as np
import pandas as pd
import Utils.MultiFLASH as multiflash
import argparse
import pathlib
import csv

BASE = pathlib.Path()

def BEfinder(user_name, project_name, sample_date, sample_list):
    sample_dir = BASE / user_name / project_name / "Input"

#ref csv file : amplicon sequence, target start position
    with open(f"{sample_dir / sample_date}.csv", 'r') as f:
        ref_file = csv.reader(f)
        ref_data = next(ref_file)

    ref_seq = ref_data[0]
    target_pos = int(ref_data[1])

    for sample in sample_list:
        f = open(f'{sample_dir / sample_date / sample / sample}.extendedFrags.fastq')
        data = f.readlines()
        seq_list = list()
        num = 0

        for line in data:
            line = line.strip()
            line = line[4:]
            if line.startswith(ref_seq[:20]):
                target=line[target_pos:target_pos+25]
                target_list = list(target)
                seq_list.append(target_list)
                num += 1


        #BEfinder ver 2.0 : result with all base 
        seq_vector = np.array(seq_list)
        ref_vector = np.array(list(ref_seq[target_pos:target_pos+25]))
        print(seq_vector.shape, ref_vector.shape)
        base_array = np.array([[['A']], [['T']], [['G']], [['C']]])
        seq_vector = seq_vector==base_array
        ref_vector = ref_vector==base_array
        mut_vector = (seq_vector!=ref_vector) & seq_vector
        print(seq_vector.shape, ref_vector.shape, mut_vector.shape)
        out = [np.array(list(ref_seq[target_pos:target_pos+25]))]
        out = np.concatenate((out, seq_vector.sum(axis=1)),axis=0)
        
        df = pd.DataFrame(out)
        df.index=[sample,'A','T','G','C']
        df.to_csv(f'{BASE / user_name / project_name / "Output" / sample_date / sample}result.csv', header=True, index=True)

if __name__ == '__main__':
    print('RUN Base Editing Finder program')

    parser = argparse.ArgumentParser(
        prog="Base Editing Finder",
        description="program for count Base Editing from fastq files",
        epilog="SKKUGE_DEV, 2023~"
    )

# command line input
    parser.add_argument(
        "-u", "--user", dest="user_name", type=str, help="The user name with no space",
    )

    parser.add_argument(
        "-p", "--project", dest="project_name", type=str, help="The project name with no space",
    )

    parser.add_argument(
        "-d", "--sample_date", dest="sample_date", type=str, help="NGS result sample_date, YYMMDD",
    )
    args = parser.parse_args()

#FLASH execute
    sample_list = multiflash.merge(args.user_name, args.project_name, args.sample_date)
    BEfinder(args.user_name, args.project_name, args.sample_date, sample_list)