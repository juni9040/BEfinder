import subprocess
import shlex
import os
import argparse
import pathlib

'''
HOW TO USE (only MultiFLASH)
1. terminal code : python FLASH_multiple_files.py -u LBJ -p test -d yymmdd
2. copy fastq file /LBJ/test/Input/yymmdd/
3. run again : python FLASH_multiple_files.py -u LBJ -p test -d yymmdd
'''

BASE = pathlib.Path()

def merge(user_name, project_name, sample_date):
#directory initialization
    input_dir = BASE / user_name / project_name / "Input"
    output_dir = BASE / user_name / project_name / "Output"
    sample_dir = input_dir / sample_date
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_dir.mkdir(parents=True, exist_ok=True)

    input_files = [file for file in os.listdir(input_dir / sample_date) if file.endswith(r'.fastq.gz')]
    input_files.sort()
    print(input_dir)

    waiting_queue = []
    for fwd, rev in zip(*[iter(input_files)] * 2):
        waiting_queue.append((fwd, rev))

    print(f"Target loaded successfully: {waiting_queue[-1]}")

    # Generate child processes for FLASH
    sample_list = []
    for fwd, rev in waiting_queue:
        filename = fwd.split("_")[1]
        sample_list.append(filename)
        subprocess.run(
            shlex.split(f"mkdir -p {input_dir / sample_date / filename}", posix=False)
        )
        cmd_input = shlex.split(
            f"flash {input_dir / sample_date / fwd} {input_dir / sample_date / rev} -M 400 -m 1 -O -o {input_dir / sample_date / filename / filename} -t 4", posix=False
        )   
        print(cmd_input)
        subprocess.run(cmd_input)
    print('flash end')

    return sample_list

if __name__ == '__main__':
    print('RUN flash multiple file program')

    parser = argparse.ArgumentParser(
        prog="FLASH_multiple_files",
        description="program for one-step FLASH fastq files",
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
    merge(args.user_name, args.project_name, args.sample_date)