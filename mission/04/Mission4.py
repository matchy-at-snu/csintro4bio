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

# Last modified: 2019-05-19

sRevCompMap = {"A": "T", "C": "G", "G": "C", "T": "A"}

from scipy.stats import fisher_exact as fisherExact
from itertools import product

sMotifList = ["".join(x) for x in product(["A", "C", "G", "T"], repeat=7)]
nMotifDict = {}
for sKey in sMotifList:
    nMotifDict[sKey] = [0, 0]  # down, notdown

global PSEUDO_COUNT
PSEUDO_COUNT = 0.0000001


class cRefSeqObj:
    def __init__(self):
        self.sGeneSym = "NULL"
        self.sNMID = "NULL"
        self.nChrID = 0  # 1~22, 23 for X, 24 for Y, 0 if not mapped correctly
        self.cStrand = "NULL"  # "+" or "-"
        self.nTxStart = 0  # transcription start position
        self.nTxEnd = 0  # transcriptin end position
        self.nCdsStart = 0  # coding region start
        self.nCdsEnd = 0  # coding region end
        self.nExonCount = 0  # number of exons
        self.nExonStarts = [0] * self.nExonCount  # exion start positions
        self.nExonEnds = [0] * self.nExonCount  # exon end positions
        self.nNumNMID = 0
        self.sExonSeq = "NULL"
        self.nExonSeqSize = 0
        self.s5UTRSeq = "NULL"
        self.n5UTRSeqSize = 0
        self.sORFSeq = "NULL"
        self.nORFSeqSize = 0
        self.s3UTRSeq = "NULL"
        self.n3UTRSeqSize = 0

    def parseRefFlatFile(self, sLine):
        sField = sLine.strip("\n").split("\t")
        sChrID = sField[2].replace("chr", "")
        sExonStarts = sField[9].rstrip(",")
        sExonEnds = sField[10].rstrip(",")
        self.sGeneSym = sField[0]
        self.sNMID = sField[1]
        self.nChrID = self.setNumChrID(sChrID)
        self.cStrand = sField[3]
        self.nTxStart = eval(sField[4])
        self.nTxEnd = eval(sField[5])
        self.nCdsStart = eval(sField[6])
        self.nCdsEnd = eval(sField[7])
        self.nExonCount = eval(sField[8])
        self.nExonStarts = [eval(x) for x in sExonStarts.split(",")]
        self.nExonEnds = [eval(x) for x in sExonEnds.split(",")]
        self.nNumNMID = eval(self.sNMID[3:].lstrip("0"))

    def getNMID(self):
        return self.sNMID

    def getChrID(self):
        return self.nChrID

    def getGeneSym(self):
        return self.sGeneSym

    def getNumNMID(self):
        return self.nNumNMID

    def getORFSeqSize(self):
        return self.nORFSeqSize

    def getORFSeq(self):
        return self.sORFSeq

    def get3UTRSeq(self):
        return self.s3UTRSeq

    def setNumChrID(self, sChrID):
        nResult = 0
        try:
            nResult = eval(sChrID)
            return nResult
        except:
            if sChrID == "X":
                nResult = 23
            elif sChrID == "Y":
                nResult = 24
            return nResult

    def setUpperGeneSym(self):
        self.sGeneSym = self.sGeneSym.upper()

    def parseTranscript(self, sChrom):
        sExonSeq = ""
        nCDSStartExon = -1
        nCDSEndExon = -1
        sORFSeq = ""
        s5UTRSeq = ""
        s3UTRSeq = ""

        for i in range(self.nExonCount):
            nStart = self.nExonStarts[i]
            nEnd = self.nExonEnds[i]
            sExonSeq += sChrom[nStart:nEnd]
            if nStart <= self.nCdsStart and nEnd >= self.nCdsStart:
                nCDSStartExon = i
            if nStart <= self.nCdsEnd and nEnd >= self.nCdsEnd:
                nCDSEndExon = i

        if nCDSStartExon != nCDSEndExon:
            for i in range(self.nExonCount):
                nStart = self.nExonStarts[i]
                nEnd = self.nExonEnds[i]
                if i < nCDSStartExon:
                    s5UTRSeq += sChrom[nStart:nEnd]
                elif i == nCDSStartExon:
                    nCdsStart = self.nCdsStart
                    s5UTRSeq += sChrom[nStart:nCdsStart]
                    sORFSeq += sChrom[nCdsStart:nEnd]
                elif i > nCDSStartExon and i < nCDSEndExon:
                    sORFSeq += sChrom[nStart:nEnd]
                elif i == nCDSEndExon:
                    nCdsEnd = self.nCdsEnd
                    sORFSeq += sChrom[nStart:nCdsEnd]
                    s3UTRSeq += sChrom[nCdsEnd:nEnd]
                elif i > nCDSEndExon:
                    s3UTRSeq += sChrom[nStart:nEnd]
        else:
            for i in range(self.nExonCount):
                nStart = self.nExonStarts[i]
                nEnd = self.nExonEnds[i]
                if i < nCDSStartExon:
                    s5UTRSeq += sChrom[nStart:nEnd]
                elif i == nCDSStartExon:
                    nCdsStart = self.nCdsStart
                    nCdsEnd = self.nCdsEnd
                    s5UTRSeq += sChrom[nStart:nCdsStart]
                    sORFSeq += sChrom[nCdsStart:nCdsEnd]
                    s3UTRSeq += sChrom[nCdsEnd:nEnd]
                elif i > nCDSEndExon:
                    s3UTRSeq += sChrom[nStart:nEnd]

        if self.cStrand == "+":
            self.sExonSeq = sExonSeq
            self.s5UTRSeq = s5UTRSeq
            self.sORFSeq = sORFSeq
            self.s3UTRSeq = s3UTRSeq
        else:
            self.sExonSeq = revComp(sExonSeq)
            self.s5UTRSeq = revComp(s3UTRSeq)
            self.sORFSeq = revComp(sORFSeq)
            self.s3UTRSeq = revComp(s5UTRSeq)

        self.nExonSeqSize = len(self.sExonSeq)
        self.nORFSeqSize = len(self.sORFSeq)
        self.n5UTRSeqSize = len(self.s5UTRSeq)
        self.n3UTRSeqSize = len(self.s3UTRSeq)


# end of class definition


def readChromSeq(sPath):
    with open(sPath, "r") as InFile:
        InFile.readline()
        sGenome = InFile.read().upper().replace("\n", "")
    return sGenome


def readRefSeq(sFileName):
    cRefSeqList = []
    with open(sFileName, "r") as InFile:
        for sLine in InFile.readlines():
            cRefSeq = cRefSeqObj()
            cRefSeq.parseRefFlatFile(sLine)
            cRefSeqList.append(cRefSeq)
    return cRefSeqList

def writeMotifList(sFileName, cMotifList):
    out_dir = pathlib.Path(__file__).parent / "solutions"
    out_dir.mkdir(exist_ok=True)
    with open(out_dir / sFileName, "w") as OutFile:
        for cMotif in cMotifList:
            OutFile.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    cMotif.getMotifSeq(),
                    cMotif.getPValue(),
                    cMotif.getMotifDown(),
                    cMotif.getNotMotifDown(),
                    cMotif.getMotifNotDown(),
                    cMotif.getNotMotifNotDown(),
                    cMotif.getRelativeRisk(),
                )
            )



def revComp(sSeq):
    sRevCompSeq = ""
    for sNuc in sSeq[::-1]:
        sRevCompSeq += sRevCompMap[sNuc]
    return sRevCompSeq


def findInternalStop(sSeq):
    sTemp = sSeq[:-3]
    for i in range(0, len(sTemp) - 2, 3):
        sCodon = sTemp[i : i + 3]
        if sCodon in ["TAA", "TGA", "TAG"]:
            return True
    return False


def removeNonNM(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        if cRefSeq.getNMID().startswith("NM_"):
            cResult.append(cRefSeq)
    return cResult


def removeAbnormalChrMap(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        if cRefSeq.getChrID() != 0:
            cResult.append(cRefSeq)
    return cResult


def removeMultiNMIDEntries(cRefSeqList):
    nCount = {}
    nDelIndex = []
    # remvoe multiple NMID
    for cRefSeq in cRefSeqList:
        sNMID = cRefSeq.getNMID()
        if sNMID not in nCount:
            nCount[sNMID] = cRefSeq
        else:
            nDelIndex.append(nCount[sNMID])
            nDelIndex.append(cRefSeq)
    nDelIndex = set(nDelIndex)
    for cRefSeq in nDelIndex:
        cRefSeqList.remove(cRefSeq)


def removeWrongLenORF(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        nORFSeqSize = cRefSeq.getORFSeqSize()
        if (nORFSeqSize % 3) == 0:
            cResult.append(cRefSeq)
    return cResult


def removeWrongCodonORF(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        sORFSeq = cRefSeq.getORFSeq()
        if not sORFSeq.startswith("ATG"):
            continue
        elif not sORFSeq[-3:] in ["TAA", "TAG", "TGA"]:
            continue
        elif findInternalStop(sORFSeq):
            continue
        else:
            cResult.append(cRefSeq)
    return cResult


def removeMultiIsoform(cRefSeqList):
    cResult = []
    nEntry = {}
    for i in range(len(cRefSeqList)):
        sGeneSym = cRefSeqList[i].sGeneSym
        if not sGeneSym in nEntry:
            nEntry[sGeneSym] = [i]
        else:
            nEntry[sGeneSym].append(i)
    for nIndexList in nEntry.values():
        nIndexList.sort(key=lambda x: cRefSeqList[x].getNumNMID())
        nIndex = nIndexList[0]
        cResult.append(cRefSeqList[nIndex])
    return cResult


class cHeLaGene:
    def __init__(self):
        self.sGeneSym = "NULL"
        self.fFoldChange = 0.0
        self.bHighlyDown = False
        self.s3UTRSeq = "NULL"

    def parseDataSet(self, sEntry):
        sField = sEntry.strip("\n").split("\t")
        self.sGeneSym = sField[0].upper()
        self.fFoldChange = eval(sField[1])
        if self.fFoldChange < -0.5:
            self.bHighlyDown = True

    def set3UTRSeq(self, s3UTRSeq):
        self.s3UTRSeq = s3UTRSeq

    def getGeneSym(self):
        return self.sGeneSym

    def getFlag(self):
        return self.bHighlyDown

    def get3UTRSeq(self):
        return self.s3UTRSeq


class cMotifObj:
    def __init__(self):
        self.sMotifSeq = "NULL"
        self.nContigTable = [[0, 0], [0, 0]]
        self.fRelativeRisk = 0.0
        self.fPValue = 0.0

    def setMotif(self, sMotifSeq):
        self.sMotifSeq = sMotifSeq

    def setContigTable(self, n1, n2, n3, n4):
        self.nContigTable[0][0] = n1
        self.nContigTable[0][1] = n2
        self.nContigTable[1][0] = n3
        self.nContigTable[1][1] = n4

    def setRelativeRisk(self):
        n1 = self.nContigTable[0][0]
        n2 = self.nContigTable[0][1]
        n3 = self.nContigTable[1][0]
        n4 = self.nContigTable[1][1]
        fA = (n1 + PSEUDO_COUNT) / (n1 + n2 + PSEUDO_COUNT)
        fB = (n3 + PSEUDO_COUNT) / (n3 + n4 + PSEUDO_COUNT)
        self.fRelativeRisk = fA / fB

    def setPValue(self):
        nTable = self.nContigTable
        fOddsRatio, fPValue = fisherExact(nTable)
        self.fPValue = fPValue

    def getMotifSeq(self):
        return self.sMotifSeq

    def getMotifDown(self):
        return self.nContigTable[0][0]

    def getNotMotifDown(self):
        return self.nContigTable[0][1]

    def getMotifNotDown(self):
        return self.nContigTable[1][0]

    def getNotMotifNotDown(self):
        return self.nContigTable[1][1]

    def getRelativeRisk(self):
        return self.fRelativeRisk

    def getPValue(self):
        return self.fPValue


def BinarySearch(arr, l, r, x):
    while l <= r:
        mid = l + (r - l) // 2
        if arr[mid].getGeneSym().upper() == x:
            return mid
        elif arr[mid].getGeneSym().upper() < x:
            l = mid + 1
        else:
            r = mid - 1
    return -1


def readDataSet(sFileName):
    cDataSet = []
    with open(sFileName, "r") as InFile:
        for sLine in InFile:
            cData = cHeLaGene()
            cData.parseDataSet(sLine)
            cDataSet.append(cData)
    return cDataSet


def removeNonExist(cDataSet):
    cResult = []
    for cGene in cDataSet:
        if cGene.get3UTRSeq() == "NULL":
            continue
        else:
            cResult.append(cGene)
    return cResult


def removeNotEnriched(cDataSet):
    cResult = []
    for cGene in cDataSet:
        if cGene.getRelativeRisk() < 1.0:
            continue
        else:
            cResult.append(cGene)
    return cResult


def main():
    # read DataSet, at the same time decide which ones are highly downregulated
    data_dir = utils.get_data_dir()

    print("Read Mission4 DataSet.")
    cDataSet = readDataSet(data_dir / "Mission4_Dataset.txt")

    print("Read ref seq.")
    sFileName = "refFlat.txt"
    sChr = "chr{}.fa"
    print("Start.")
    cRefSeqList = readRefSeq(data_dir / sFileName)

    # Task 1. Read results into cRefSeqList, count the number of entries in the RefFalt file.
    print("The number of entries: ", len(cRefSeqList))

    # Task 2. Remove non-NM sequences and keep the ones aligned only on chr 1-22, X, or Y. Count the number of entries left.
    print(
        (
            "Remove non-NM sequences "
            "and keep the ones aligned only on chr1-22, X or Y."
        )
    )
    cRefSeqList = removeNonNM(cRefSeqList)
    print("The number of entries left after removing non-NM: ", len(cRefSeqList))
    cRefSeqList = removeAbnormalChrMap(cRefSeqList)
    print(
        "The number of entries left after removing wrong mappings: ", len(cRefSeqList)
    )

    # Task 3. Remove NM sequences that have multiple entries in the RefFlat file. Count the number of entries left.
    print("Remove sequences with multiple entries.")
    removeMultiNMIDEntries(cRefSeqList)
    print("The number of entries left after removing multiple NMID:", len(cRefSeqList))

    # Task 4. Remove NM sequences that have wrong ORFs: no start codon, no stop codon, ORF size of 3N+1 or 3N+2, or internal stop codons. Count the number of entries left
    path = utils.get_species_genome_dir("human")
    for i in range(1, 25):
        if i == 23:
            sChrom = readChromSeq(path / sChr.format("X"))
        elif i == 24:
            sChrom = readChromSeq(path / sChr.format("Y"))
        else:
            sChrom = readChromSeq(path / sChr.format(i))
        print("Start to parse transcriptome for chrom {}".format(i))
        for cRefSeq in cRefSeqList:
            if cRefSeq.nChrID == i:
                cRefSeq.parseTranscript(sChrom)
    del sChrom

    print("Remove sequences with ORFs of wrong lengths.")
    cRefSeqList = removeWrongLenORF(cRefSeqList)
    print("The number of enries left after removing wrong ORF:", len(cRefSeqList))
    print("Remove sequences with ORFs of wrong codon.")
    cRefSeqList = removeWrongCodonORF(cRefSeqList)
    print("The number of enries left after removing wrong ORF:", len(cRefSeqList))

    # Task 5. For each gene that has multiple isoforms, select a representative isoform by picking up a RefSeq sequence that has the lowest NM ID(an integer). Count the number of entries left.
    print("Reduce multiple isoforms.")
    cRefSeqList = removeMultiIsoform(cRefSeqList)
    print(
        "The number of enries left after removing multiple isoform:", len(cRefSeqList)
    )

    # # Task 6. Sort the list of the RefSeq class objects by the ascending order of the RefSeq ID (not a string but an integer), and store the RefSeq ID and the Gene Symbol in the provided template file.
    for cRefSeq in cRefSeqList:
        cRefSeq.setUpperGeneSym()

    cRefSeqList.sort(key=lambda x: x.getGeneSym())

    ############################# Start of Mission 4 #############################

    # Add 3UTR Seq to
    print("Start to add 3UTR to HeLaGene")

    for cGene in cDataSet:
        x = cGene.getGeneSym()
        i = BinarySearch(cRefSeqList, 0, len(cRefSeqList) - 1, x)
        if i < 0:
            continue
        else:
            cGene.set3UTRSeq(cRefSeqList[i].get3UTRSeq())

    cDataSet = removeNonExist(cDataSet)

    # Built Contingency Table for each motif

    nTotalDown = 0
    nTotalNotDown = 0
    print("Building contigency table...")
    for cGene in cDataSet:
        s3UTRSeq = cGene.get3UTRSeq()
        nLength = len(s3UTRSeq)
        if nLength < 7:
            if cGene.getFlag():
                nTotalDown += 1
            else:
                nTotalNotDown += 1
            continue
        sMotifSet = []
        for i in range(nLength - 6):
            sMotif = s3UTRSeq[i : i + 7]
            sMotifSet.append(sMotif)
        sMotifSet = set(sMotifSet)
        if cGene.getFlag():
            nTotalDown += 1
            for x in sMotifSet:
                nMotifDict[x][0] += 1
        else:
            nTotalNotDown += 1
            for x in sMotifSet:
                nMotifDict[x][1] += 1

    cMotifList = []
    print("Doing Fisher Exact Test...")
    for sKey in nMotifDict:
        nMotifDown = nMotifDict[sKey][0]
        nMotifNotDown = nMotifDict[sKey][1]
        cMotif = cMotifObj()
        cMotif.setMotif(sKey)
        cMotif.setContigTable(
            nMotifDown,
            nTotalDown - nMotifDown,
            nMotifNotDown,
            nTotalNotDown - nMotifNotDown,
        )
        cMotif.setRelativeRisk()
        cMotif.setPValue()
        cMotifList.append(cMotif)

    cMotifList_Enriched = removeNotEnriched(cMotifList)

    cMotifList_Enriched.sort(key=lambda x: x.getPValue())

    # Write only top 10
    writeMotifList("Mission4_MotifList.txt", cMotifList_Enriched[:10])


if __name__ == "__main__":
    main()
