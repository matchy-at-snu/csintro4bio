import os
import sys
import pathlib

# Get the path to the mission directory
mission_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Add the mission directory to the Python path
sys.path.insert(0, mission_dir)

# Import the path modifier
import path_modifier
from mission import utils

# Last modified: 2019-04-14

# A program to perform analysis on 22 autosomes for human, 29 chromosomes(1~28, 32) for chicken, 5(2R, 2L, 3R, 3L & 4) for fruit fly, 5(simply 1~5) for C. elegans


def readSeq(sPath):
    with open(sPath, "r") as InFile:
        InFile.readline()
        sGenome = InFile.read().upper().replace("\n", "")
    return sGenome
    # end of readSeq


def writeResult(
    sSpecies, nMonoNtDict, fMonoFreqDict, nDiNtDict, fExpDiFreqDict, fDiFreqDict
):
    # output the results into ./solutions/
    out_dir = pathlib.Path(__file__) / "solutions"
    out_dir.mkdir(exist_ok=True)
    sFileName = "{}genome.solutions.txt".format(sSpecies)

    with open(out_dir / sFileName, "w") as OutFile:
        sMonoNtKeyList = sorted(list(nMonoNtDict.keys()))
        for sNuc in sMonoNtKeyList:
            nNum = nMonoNtDict[sNuc]
            fFreq = fMonoFreqDict[sNuc]
            sWriteString = ("{}: num = {},\n" "   freq = {};\n").format(
                sNuc, nNum, fFreq
            )
            OutFile.write(sWriteString)
        OutFile.write(
            "-----------------------------------------------------------------------------\n"
        )
        sDiNtKeyList = sorted(list(nDiNtDict.keys()))
        for sSeq in sDiNtKeyList:
            nNum = nDiNtDict[sSeq]
            fExpFreq = fExpDiFreqDict[sSeq]
            fFreq = fDiFreqDict[sSeq]
            sWriteString = (
                "{}: num = {},\n" "    expected freq = {},\n" "    freq = {};\n"
            ).format(sSeq, nNum, fExpFreq, fFreq)
            OutFile.write(sWriteString)
    print("Please check {} for the output results.".format(sFileName))
    # end of writeResult


def monoNtCount(nResult, sSeq):
    for sNuc in sSeq:
        if not sNuc in nResult:
            nResult[sNuc] = 1
        else:
            nResult[sNuc] += 1
    # end of monoNtCount


def diNtCount(nResult, sSeq):
    n = len(sSeq)
    for i in range(0, n - 1):
        sSubSeq = sSeq[i : i + 2]
        if not sSubSeq in nResult:
            nResult[sSubSeq] = 1
        else:
            nResult[sSubSeq] += 1
    # end of diNtCount


def withoutN(nNtDict):
    sKeyList = list(nNtDict.keys())
    for sKey in sKeyList:
        if "N" in sKey:
            del nNtDict[sKey]
    # end of withoutN


def freqCalc(nNtDict):
    fResult = {}
    nTotal = sum(nNtDict.values())
    for sNuc in nNtDict:
        fResult[sNuc] = nNtDict[sNuc] / nTotal
    return fResult
    # end of freqCalc


def expFreqCalc(fMonoNtFreqDict, nDiNtDict):
    fResult = {}
    for sSeq in nDiNtDict:
        fResult[sSeq] = fMonoNtFreqDict[sSeq[0]] * fMonoNtFreqDict[sSeq[1]]
    return fResult
    # end of expFreqCalc


def pathRedirect(sSpecies):
    sPath = utils.get_species_genome_dir(sSpecies)
    if sSpecies == "human":
        rgChrList = list(range(1, 23))
    elif sSpecies == "chicken":
        rgChrList = list(range(1, 29)) + [32]
    elif sSpecies == "fruitfly":
        rgChrList = ["2R", "2L", "3R", "3L", "4"]
    else:
        rgChrList = ["I", "II", "III", "IV", "V"]
    return sPath, rgChrList
    # end of pathRedirect


def chrAnalysis(sSpecies):
    # 1. path initialization according to input species
    sPath, rgChrList = pathRedirect(sSpecies)
    # 2. count the number of mononucleotides and the number of di-nucleotides
    nMonoNtDict = {}
    nDiNtDict = {}
    for i in rgChrList:
        print("Analysing {} genome chrom {}".format(sSpecies, i))
        sChrom = readSeq(sPath / f"chr{i}.fa")
        monoNtCount(nMonoNtDict, sChrom)
        print("Mononucleotides analysis finished.")
        diNtCount(nDiNtDict, sChrom)
        print("Dinucleotides analysis finished.")
    # 2.1 Get rid of the ambiguous Ns
    print("Getting rid of the ambiguous Ns...")
    withoutN(nMonoNtDict)
    withoutN(nDiNtDict)
    # 3. calculate the relative frequencies of mononucleotides
    print("Calculating mononucleotides relative frequencies...")
    fMonoFreqDict = freqCalc(nMonoNtDict)
    # 4. calculate the expected relative frequencies of dinucleotides using the frequencies of coresponding mononucleotides
    print("Calculating dinucleotides expected frequencies...")
    fExpDiFreqDict = expFreqCalc(fMonoFreqDict, nDiNtDict)
    # 5. calculate the real frequencies of dinucleotides
    print("Calculating dinucleotides real frequencies...")
    fDiFreqDict = freqCalc(nDiNtDict)
    # 6. write the results into OutFile
    print("Writing results...")
    writeResult(
        sSpecies, nMonoNtDict, fMonoFreqDict, nDiNtDict, fExpDiFreqDict, fDiFreqDict
    )
    # end of chrAnalysis


def main():
    while True:
        # Get Input
        sSpecies = input(
            (
                "Please enter the species name that you want to analyse:\n"
                "Species available:\n"
                "human, chicken, Drosophila melanogaster(enter 'fruitfly'), C.elegans(enter 'celegans')\n"
                "Enter 'end' if you would like to end the analysis.\n"
                "The results will be stored in txt files.\n"
            )
        )
        # Process
        if sSpecies == "end":
            break
        chrAnalysis(sSpecies)
    print("Genome analysis finished.")
    # end of main


if __name__ == "__main__":
    main()
