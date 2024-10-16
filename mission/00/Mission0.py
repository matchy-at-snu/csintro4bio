import os
import sys

# Get the path to the mission directory
mission_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add the mission directory to the Python path
sys.path.insert(0, mission_dir)

# Import the path modifier
import path_modifier
from mission import utils

# Your mission code goes here


from mission.utils import *

# Last modified: 2019-04-11


def ntCount(result, seq):
    for nt in seq:
        if not nt in result:
            continue
        result[nt] += 1


def main():
    # read chr2.fa file
    with open(utils.get_data_dir() / "hg38" / "chroms"/ "chr2.fa") as file:
        file.readline()
        genome = file.read().upper().replace("\n", "")

    total = len(genome)
    genome_mod = genome.replace("N", "")
    total_nonHetero = len(genome_mod)

    print("The total length of the genome is {}".format(total))
    print(
        "The total length of the genome without the unsequenced heterochromatin is {}".format(
            total_nonHetero
        )
    )

    print("Result:")
    ntDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    ntCount(ntDict, genome_mod)
    for key in ntDict:
        print("{}: {}".format(key, ntDict[key]))
        freq = ntDict[key] / total_nonHetero * 100
        print("The frequency of {} is: {} %".format(key, freq))


if __name__ == "__main__":
    main()
