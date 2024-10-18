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

# Last modified: 2019-05-29

################################## Mission 5 ##################################

# ============================= Import Package ============================== #

from scipy.stats import fisher_exact as fisherExact
from itertools import product
from statsmodels.stats.multitest import multipletests as mt

# =========================================================================== #

# ============================ Global Variables ============================= #

sRevCompMap = {"A": "T", "C": "G", "G": "C", "T": "A"}
sRNARevCompMap = {"A": "U", "C": "G", "G": "C", "T": "A"}
sMotifList = ["".join(x) for x in product(["A", "C", "G", "T"], repeat=7)]

global PSEUDO_COUNT
PSEUDO_COUNT = 0.0000001

# =========================================================================== #

# ============================ Class Definition ============================= #


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
        self.nExonStarts = [0] * self.nExonCount  # exon start positions
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

    def get5UTRSeq(self):
        return self.s5UTRSeq

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


### end of class definition


class cHeLaGene:
    def __init__(self):
        self.sGeneSym = "NULL"
        self.fFoldChange = 0.0
        self.bHighlyDown = False
        self.sORFSeq = "NULL"
        self.s3UTRSeq = "NULL"
        self.s5UTRSeq = "NULL"

    def parseDataSet(self, sEntry):
        sField = sEntry.strip("\n").split("\t")
        self.sGeneSym = sField[0].upper()
        self.fFoldChange = eval(sField[1])
        if self.fFoldChange < -0.5:
            self.bHighlyDown = True

    def setORFSeq(self, sORFSeq):
        self.sORFSeq = sORFSeq

    def set3UTRSeq(self, s3UTRSeq):
        self.s3UTRSeq = s3UTRSeq

    def set5UTRSeq(self, s5UTRSeq):
        self.s5UTRSeq = s5UTRSeq

    def getGeneSym(self):
        return self.sGeneSym

    def getFlag(self):
        return self.bHighlyDown

    def getORFSeq(self):
        return self.sORFSeq

    def get3UTRSeq(self):
        return self.s3UTRSeq

    def get5UTRSeq(self):
        return self.s5UTRSeq


### end of class definition


class cMotifObj:
    def __init__(self):
        self.sMotifSeq = "NULL"
        self.nContigTable = [[0, 0], [0, 0]]
        self.fRelativeRisk = 0.0
        self.fPValue = 0.0
        self.sMiRType = ""
        self.sMiRMatch = []

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

    def correctPVal(self, fPValue):
        self.fPValue = fPValue

    def setMiRMatch(self, sName, sList):
        self.sMiRType = sName
        self.sMiRMatch += sList

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

    def getMiRType(self):
        return self.sMiRType

    def getMiRMatch(self):
        return ",".join(self.sMiRMatch)


### end of class definition


class cMiRObj:
    def __init__(self):
        sMiRName = "NULL"
        sMiRSeq = "NULL"

    def setMiRName(self, sMiRName):
        self.sMiRName = sMiRName

    def setMiRSeq(self, sMiRSeq):
        self.sMiRSeq = sMiRSeq

    def getMiRName(self):
        return self.sMiRName

    def getMiRSeq(self):
        return self.sMiRSeq


# =========================================================================== #

# =========================== Function Definition =========================== #
# ---------------------------- Useful Alogrithms ---------------------------- #


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


### end of Binary Search

# --------------------------------------------------------------------------- #


def readChromSeq(sPath):
    with open(sPath, "r") as InFile:
        InFile.readline()
        sGenome = InFile.read().upper().replace("\n", "")
    return sGenome


### end of function definition


def readRefSeq(sFileName):
    cRefSeqList = []
    with open(sFileName, "r") as InFile:
        for sLine in InFile.readlines():
            cRefSeq = cRefSeqObj()
            cRefSeq.parseRefFlatFile(sLine)
            cRefSeqList.append(cRefSeq)
    return cRefSeqList


### end of function definition


def readDataSet(sFileName):
    cDataSet = []
    with open(sFileName, "r") as InFile:
        for sLine in InFile:
            cData = cHeLaGene()
            cData.parseDataSet(sLine)
            cDataSet.append(cData)
    return cDataSet


### end of function definition


def readMiRDB(sFileName):
    cMiRDataSet = []
    with open(sFileName, "r") as InFile:
        sSpecies = "NULL"
        sMiRName = "NULL"
        for sLine in InFile:
            if sLine.startswith(">"):
                sField = sLine.rstrip("\n").split(" ")
                sSpecies = sField[2] + " " + sField[3]
                sMiRName = sField[4]
            else:
                if sSpecies != "Homo sapiens":
                    continue
                else:
                    cMiRNA = cMiRObj()
                    cMiRNA.setMiRName(sMiRName)
                    cMiRNA.setMiRSeq(sLine.rstrip("\n"))
                    cMiRDataSet.append(cMiRNA)
    return cMiRDataSet


def revComp(sSeq, bRNA=False):
    sRevCompSeq = ""
    if bRNA:
        for sNuc in sSeq[::-1]:
            sRevCompSeq += sRNARevCompMap[sNuc]
    else:
        for sNuc in sSeq[::-1]:
            sRevCompSeq += sRevCompMap[sNuc]
    return sRevCompSeq


### end of function definition


def findInternalStop(sSeq):
    sTemp = sSeq[:-3]
    for i in range(0, len(sTemp) - 2, 3):
        sCodon = sTemp[i : i + 3]
        if sCodon in ["TAA", "TGA", "TAG"]:
            return True
    return False


### end of function definition


def removeNonNM(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        if cRefSeq.getNMID().startswith("NM_"):
            cResult.append(cRefSeq)
    return cResult


### end of function definition


def removeAbnormalChrMap(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        if cRefSeq.getChrID() != 0:
            cResult.append(cRefSeq)
    return cResult


### end of function definition


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


### end of function definition


def removeWrongLenORF(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        nORFSeqSize = cRefSeq.getORFSeqSize()
        if (nORFSeqSize % 3) == 0:
            cResult.append(cRefSeq)
    return cResult


### end of function definition


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


### end of function definition


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


### end of function definition


def remove3UTRNonExist(cDataSet):
    cResult = []
    for cGene in cDataSet:
        if cGene.get3UTRSeq() == "NULL":
            continue
        else:
            cResult.append(cGene)
    return cResult


### end of function definition


def removeORFNonExist(cDataSet):
    cResult = []
    for cGene in cDataSet:
        if cGene.getORFSeq() == "NULL":
            continue
        else:
            cResult.append(cGene)
    return cResult


### end of function definition


def remove5UTRNonExist(cDataSet):
    cResult = []
    for cGene in cDataSet:
        if cGene.get5UTRSeq() == "NULL":
            continue
        else:
            cResult.append(cGene)
    return cResult


### end of function definition


def removeNotEnriched(cDataSet):
    cResult = []
    for cGene in cDataSet:
        if cGene.getRelativeRisk() < 1.0:
            continue
        else:
            cResult.append(cGene)
    return cResult


### end of function definition


def pValCorrection(cMotifList):
    fPVal = [x.getPValue() for x in cMotifList]
    bReject, fCPVal, alphacSidak, alphacBonf = mt(fPVal, method="bonferroni")
    for i in range(len(cMotifList)):
        cMotifList[i].correctPVal(fCPVal[i])


### end of function definition


def prepRefSeqList(sFileName):
    print("Read ref seq.")
    sFileName = sFileName
    sPath = "chr{}.fa"
    print("Start.")
    cRefSeqList = readRefSeq(sFileName)
    print(
        (
            "Remove non-NM sequences "
            "and keep the ones aligned only on chr1-22, X or Y."
        )
    )
    cRefSeqList = removeNonNM(cRefSeqList)
    cRefSeqList = removeAbnormalChrMap(cRefSeqList)
    print("Remove sequences with multiple entries.")
    removeMultiNMIDEntries(cRefSeqList)

    path = utils.get_species_genome_dir("human")
    for i in range(1, 25):
        if i == 23:
            sChrom = readChromSeq(path / sPath.format("X"))
        elif i == 24:
            sChrom = readChromSeq(path / sPath.format("Y"))
        else:
            sChrom = readChromSeq(path / sPath.format(i))
        print("Start to parse transcriptome for chrom {}".format(i))
        for cRefSeq in cRefSeqList:
            if cRefSeq.nChrID == i:
                cRefSeq.parseTranscript(sChrom)
    del sChrom

    print("Remove sequences with ORFs of wrong lengths.")
    cRefSeqList = removeWrongLenORF(cRefSeqList)
    print("Remove sequences with ORFs of wrong codon.")
    cRefSeqList = removeWrongCodonORF(cRefSeqList)
    print("Reduce multiple isoforms.")
    cRefSeqList = removeMultiIsoform(cRefSeqList)

    # Make every gene symbol to upper case
    for cRefSeq in cRefSeqList:
        cRefSeq.setUpperGeneSym()

    # Sort the cRefSeqList by the order of GeneSymbol (ASCII, ascending manner)
    cRefSeqList.sort(key=lambda x: x.getGeneSym())
    return cRefSeqList


### end of function definition


def do3UTRAnalysis(cRefSeqList, cDataSet):
    nMotifDict = {}
    for sKey in sMotifList:
        nMotifDict[sKey] = [0, 0]  # down, notdown
    print("Preserve only entries occured in both refFlat file and dataset.")
    for cGene in cDataSet:
        x = cGene.getGeneSym()
        i = BinarySearch(cRefSeqList, 0, len(cRefSeqList) - 1, x)
        if i < 0:
            continue
        else:
            cGene.set3UTRSeq(cRefSeqList[i].get3UTRSeq())
    cDataSet = remove3UTRNonExist(cDataSet)

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

    c3UTRList = []
    print("Doing Fisher Exact Test...")
    for sKey in nMotifDict:
        nMotifDown = nMotifDict[sKey][0]
        nMotifNotDown = nMotifDict[sKey][1]
        c3UTR = cMotifObj()
        c3UTR.setMotif(sKey)
        c3UTR.setContigTable(
            nMotifDown,
            nTotalDown - nMotifDown,
            nMotifNotDown,
            nTotalNotDown - nMotifNotDown,
        )
        c3UTR.setRelativeRisk()
        c3UTR.setPValue()
        c3UTRList.append(c3UTR)

    c3UTRList.sort(key=lambda x: x.getPValue())
    pValCorrection(c3UTRList)
    c3UTRList = removeNotEnriched(c3UTRList)

    cResult = []
    for i in range(5):
        cObj = c3UTRList[i]
        cResult.append(cObj)
    del cDataSet
    del c3UTRList
    return cResult


### end of function definition


def doORFAnalysis(cRefSeqList, cDataSet):
    nMotifDict = {}
    for sKey in sMotifList:
        nMotifDict[sKey] = [0, 0]  # down, notdown
    for cGene in cDataSet:
        x = cGene.getGeneSym()
        i = BinarySearch(cRefSeqList, 0, len(cRefSeqList) - 1, x)
        if i < 0:
            continue
        else:
            cGene.setORFSeq(cRefSeqList[i].getORFSeq())
    cDataSet = removeORFNonExist(cDataSet)

    # Built Contingency Table for each motif
    nTotalDown = 0
    nTotalNotDown = 0
    print("Building contigency table...")
    for cGene in cDataSet:
        sORFSeq = cGene.getORFSeq()
        nLength = len(sORFSeq)
        if nLength < 7:
            if cGene.getFlag():
                nTotalDown += 1
            else:
                nTotalNotDown += 1
            continue
        sMotifSet = []
        for i in range(nLength - 6):
            sMotif = sORFSeq[i : i + 7]
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

    cORFList = []
    print("Doing Fisher Exact Test...")
    for sKey in nMotifDict:
        nMotifDown = nMotifDict[sKey][0]
        nMotifNotDown = nMotifDict[sKey][1]
        cORF = cMotifObj()
        cORF.setMotif(sKey)
        cORF.setContigTable(
            nMotifDown,
            nTotalDown - nMotifDown,
            nMotifNotDown,
            nTotalNotDown - nMotifNotDown,
        )
        cORF.setRelativeRisk()
        cORF.setPValue()
        cORFList.append(cORF)

    cORFList.sort(key=lambda x: x.getPValue())
    pValCorrection(cORFList)
    cORFList = removeNotEnriched(cORFList)
    cResult = []
    for i in range(5):
        cObj = cORFList[i]
        cResult.append(cObj)
    del cDataSet
    del cORFList
    return cResult


### end of function definition


def do5UTRAnalysis(cRefSeqList, cDataSet):
    nMotifDict = {}
    for sKey in sMotifList:
        nMotifDict[sKey] = [0, 0]  # down, notdown
    print("Preserve only entries occured in both refFlat file and dataset.")
    for cGene in cDataSet:
        x = cGene.getGeneSym()
        i = BinarySearch(cRefSeqList, 0, len(cRefSeqList) - 1, x)
        if i < 0:
            continue
        else:
            cGene.set5UTRSeq(cRefSeqList[i].get5UTRSeq())
    cDataSet = remove5UTRNonExist(cDataSet)

    # Built Contingency Table for each motif
    nTotalDown = 0
    nTotalNotDown = 0
    print("Building contigency table...")
    for cGene in cDataSet:
        s5UTRSeq = cGene.get5UTRSeq()
        nLength = len(s5UTRSeq)
        if nLength < 7:
            if cGene.getFlag():
                nTotalDown += 1
            else:
                nTotalNotDown += 1
            continue
        sMotifSet = []
        for i in range(nLength - 6):
            sMotif = s5UTRSeq[i : i + 7]
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

    c5UTRList = []
    print("Doing Fisher Exact Test...")
    for sKey in nMotifDict:
        nMotifDown = nMotifDict[sKey][0]
        nMotifNotDown = nMotifDict[sKey][1]
        c5UTR = cMotifObj()
        c5UTR.setMotif(sKey)
        c5UTR.setContigTable(
            nMotifDown,
            nTotalDown - nMotifDown,
            nMotifNotDown,
            nTotalNotDown - nMotifNotDown,
        )
        c5UTR.setRelativeRisk()
        c5UTR.setPValue()
        c5UTRList.append(c5UTR)

    c5UTRList.sort(key=lambda x: x.getPValue())
    pValCorrection(c5UTRList)
    c5UTRList = removeNotEnriched(c5UTRList)

    cResult = []
    for i in range(5):
        cObj = c5UTRList[i]
        cResult.append(cObj)
    del cDataSet
    del c5UTRList
    return cResult


### end of function definition


def miRType(sMotif, cMiRDataSet):
    # 7mer-A1 branch
    if sMotif.endswith("A"):
        sName = "7mer-A1"
        sList = []
        sSeed = revComp(sMotif, True)
        for cMiR in cMiRDataSet:
            sMiRNA = cMiR.getMiRSeq()
            sMiRName = cMiR.getMiRName()
            if sMiRNA[1:7] == sSeed[1:]:
                sList.append(sMiRName)
        if len(sList) != 0:
            return sName, sList
        else:
            return "", ""
    # 7mer-m8 branch
    else:
        sName = "7mer-m8"
        sList = []
        sSeed = revComp(sMotif, True)
        for cMiR in cMiRDataSet:
            sMiRNA = cMiR.getMiRSeq()
            sMiRName = cMiR.getMiRName()
            if sMiRNA[1:8] == sSeed:
                sList.append(sMiRName)
        if len(sList) != 0:
            return sName, sList
        else:
            return "", ""


### end of function definition

# =========================================================================== #

# ============================== Main Function ============================== #


def main():
    print("Prepare the RefSeqList.")
    data_dir = utils.get_data_dir()
    cRefSeqList = prepRefSeqList(data_dir / "refFlat.txt")
    sFileName = "Mission5_Dataset{}.txt"

    cDataSetList = []
    for i in range(1, 4):
        cDataSet = readDataSet(data_dir / sFileName.format(i))
        cDataSetList.append(cDataSet)

    cMiRDataSet = readMiRDB(data_dir / "mature.fa")

    cResultTuple = []
    for cDataSet in cDataSetList:
        print("Doing 3UTR Analysis")
        c3UTR = do3UTRAnalysis(cRefSeqList, cDataSet)
        print("Doing ORF Analysis")
        cORF = doORFAnalysis(cRefSeqList, cDataSet)
        cResultTuple.append((c3UTR, cORF))

    print("Matching miRNA type")
    num = 1
    for cTuple in cResultTuple:
        print("\n")
        print("The 3UTR results for Mission5_Dataset_{}.txt: ".format(num))
        out_dir = pathlib.Path(__file__).parent / "solutions"
        out_dir.mkdir(exist_ok=True)
        outfname = "mission5_dataset{}_{}.result.txt"
        c3UTR = cTuple[0]
        with open(out_dir / outfname.format(num, "3utr"), "w") as f:
            for cMotif in c3UTR:
                sName, sMatchList = miRType(cMotif.getMotifSeq(), cMiRDataSet)
                cMotif.setMiRMatch(sName, sMatchList)
                sLine = (
                    "\t".join(
                        map(
                            str,
                            [
                                cMotif.getMotifSeq(),
                                cMotif.getPValue(),
                                cMotif.getMotifDown(),
                                cMotif.getNotMotifDown(),
                                cMotif.getMotifNotDown(),
                                cMotif.getNotMotifNotDown(),
                                cMotif.getRelativeRisk(),
                                cMotif.getMiRType() + " " + cMotif.getMiRMatch(),
                            ],
                        )
                    )
                    + "\n"
                )
                f.write(sLine)
                print(sLine, end="")
        print("\n")
        print("The ORF results for Mission5_Dataset_{}.txt: ".format(num))
        cORF = cTuple[1]
        with open(out_dir / outfname.format(num, "orf"), "w") as f:
            for cMotif in cORF:
                sName, sMatchList = miRType(cMotif.getMotifSeq(), cMiRDataSet)
                cMotif.setMiRMatch(sName, sMatchList)
                sLine = (
                    "\t".join(
                        map(
                            str,
                            [
                                cMotif.getMotifSeq(),
                                cMotif.getPValue(),
                                cMotif.getMotifDown(),
                                cMotif.getNotMotifDown(),
                                cMotif.getMotifNotDown(),
                                cMotif.getNotMotifNotDown(),
                                cMotif.getRelativeRisk(),
                                cMotif.getMiRType() + " " + cMotif.getMiRMatch(),
                            ],
                        )
                    )
                    + "\n"
                )
                f.write(sLine)
                print(sLine, end="")
        num += 1

    # ADDITIONAL ANALYSIS ON 5UTR

    # cResult = []
    # for cDataSet in cDataSetList:
    #     print("Doing 5UTR Analysis")
    #     c5UTR = do5UTRAnalysis(cRefSeqList, cDataSet)
    #     cResult.append(c5UTR)

    # print("Matching miRNA type")
    # num = 1
    # for c5UTR in cResult:
    #     print("\n")
    #     print("The 5UTR results for Mission5_Dataset_{}.txt: ".format(num))
    #     for cMotif in c5UTR:
    #         sName, sMatchList = miRType(cMotif.getMotifSeq(), cMiRDataSet)
    #         cMotif.setMiRMatch(sName, sMatchList)
    #         print(cMotif.getMotifSeq(),     cMotif.getPValue(),
    #               cMotif.getMotifDown(),    cMotif.getNotMotifDown(),
    #               cMotif.getMotifNotDown(), cMotif.getNotMotifNotDown(),
    #               cMotif.getRelativeRisk(),
    #               cMotif.getMiRType()+" "+cMotif.getMiRMatch(),
    #               sep = "\t", end = "\n")
    #     print("\n")
    #     num += 1

    # 7mer-m8: exact match starts from the 8th nt
    # 7mer-A1: A1@3' on mRNA (Y1@miRNA, NO NEED TO MATCH)


# =========================================================================== #

if __name__ == "__main__":
    main()
