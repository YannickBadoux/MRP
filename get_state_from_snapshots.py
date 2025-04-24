from amuse.io import read_set_from_file

from analyse_result import bound

import glob
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get the states of folder of snapshots of SPM scattering experiments.")
    parser.add_argument('--path', type=str, help='Path to the folder containing the snapshots')
    parser.add_argument('--output', type=str, help='Path to the output file')

    args = parser.parse_args()
    path = args.path
    output = args.output

    # Get the list of all files in the directory
    files = glob.glob(path+'/*.amuse')