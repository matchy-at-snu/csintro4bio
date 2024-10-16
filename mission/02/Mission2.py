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

# Last modified: 2019-04-21

class cRefSeqObj:
    def __init__(self):
        self.sGeneSym =    "NULL"
        self.sNMID =       "NULL"
        self.nChrID =      0               # 1~22, 23 for X, 24 for Y, 0 if not mapped correctly
        self.cStrand =     "NULL"          # "+" or "-"
        self.nTxStart =    0               # transcription start position
        self.nTxEnd =      0               # transcriptin end position
        self.nCdsStart =   0               # coding region start
        self.nCdsEnd =     0               # coding region end
        self.nExonCount =  0               # number of exons
        self.nExonStarts = [0] * self.nExonCount # exion start positions
        self.nExonEnds =   [0] * self.nExonCount # exon end positions
        self.nNumNMID =    0

    def parseRefFlatFile(self, sLine):
        sField = sLine.strip("\n").split("\t")
        sChrID = sField[2].replace("chr","")
        sExonStarts = sField[9].rstrip(",")
        sExonEnds = sField[10].rstrip(",")
        self.sGeneSym =    sField[0]
        self.sNMID =       sField[1]
        self.nChrID =      self.getNumChrID(sChrID)
        self.cStrand =     sField[3]
        self.nTxStart =    eval(sField[4])
        self.nTxEnd =      eval(sField[5])
        self.nCdsStart =   eval(sField[6])
        self.nCdsEnd =     eval(sField[7])
        self.nExonCount =  eval(sField[8])
        self.nExonStarts = [ eval(x) for x in sExonStarts.split(",") ]
        self.nExonEnds =   [ eval(x) for x in sExonEnds.split(",") ]
        self.nNumNMID =    eval(self.sNMID[3:].lstrip("0"))

    def getNMID(self):
        return self.sNMID

    def getChrID(self):
        return self.nChrID

    def getGeneSym(self):
        return self.sGeneSym

    def getNumNMID(self):
        return self.nNumNMID

    def getNumChrID(self, sChrID):
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
    # end of class definition

def readFile(sFileName):
    cRefSeqList = []
    with open(sFileName, "r") as InFile:
        for sLine in InFile.readlines():
            cRefSeq = cRefSeqObj()
            cRefSeq.parseRefFlatFile(sLine)
            cRefSeqList.append(cRefSeq)
    return cRefSeqList

def writeRefSeq(cRefSeqList, sFileName = "mission2.result.txt"):
    out_dir = pathlib.Path(__file__).parent / "solutions"
    out_dir.mkdir(exist_ok = True)
    with open(out_dir / sFileName,"w") as OutFile:
        for cRefSeq in cRefSeqList:
            print(cRefSeq.sNMID, cRefSeq.sGeneSym,
                  sep = "\t", end = "\n", file = OutFile)


def removeNonNM(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        if cRefSeq.getNMID().startswith("NM_"):
            cResult.append(cRefSeq)
    return cResult

def removeAbnormalChrMap(cRefSeqList):
    cResult = []
    for cRefSeq in cRefSeqList:
        if cRefSeq.getChrID() != 0 :
            cResult.append(cRefSeq)
    return cResult

def removeMultiGeneSymEntries(cRefSeqList):
    cCount = {}
    nDelIndex = []
    for cRefSeq in cRefSeqList:
        sGeneSym = cRefSeq.getGeneSym()
        if sGeneSym not in cCount:
            cCount[sGeneSym] = cRefSeq
        else:
            nDelIndex.append(cCount[sGeneSym])
            nDelIndex.append(cRefSeq)
    nDelIndex = set(nDelIndex)
    for cRefSeq in nDelIndex:
        cRefSeqList.remove(cRefSeq)

def main():
    input_dir = utils.get_data_dir()
    sFileName = "refFlat.txt"
    print("Start.")
    # Task 1. Read results into cRefSeqList, count the number of entries in the RefFalt file.
    cRefSeqList = readFile(input_dir / sFileName)
    print("The number of entries: ", len(cRefSeqList))

    # Task 2. Remove non-NM sequences and keep the ones aligned only on chr 1-22, X, or Y. Count the number of entries left.
    print(("Remove non-NM sequences "
           "and keep the ones aligned only on chr1-22, X or Y."))
    cRefSeqList = removeNonNM(cRefSeqList)
    print("The number of entries left after removing non-NM: ",
          len(cRefSeqList))
    cRefSeqList = removeAbnormalChrMap(cRefSeqList)
    print("The number of entries left after removing wrong mappings: ",
          len(cRefSeqList))

    # Task 3. Remove NM sequences that have multiple entries in the RefFlat file. Count the number of entries left.
    print("Remove sequences with multiple entries.")
    removeMultiGeneSymEntries(cRefSeqList)
    print("The number of entries left after removing multiple gene symbol:",
           len(cRefSeqList))

    # Task 4. Sort the list of the RefSeq class objects by the ascending order of the RefSeq ID (not a string but an integer), and store the RefSeq ID and the Gene Symbol in the provided template file.
    cRefSeqList.sort(key = lambda x: x.getNumNMID())
    writeRefSeq(cRefSeqList)


if __name__ == "__main__":
    main()
