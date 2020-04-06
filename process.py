#!/usr/bin/env python3
import atexit
import glob
import os
import re
import sys
import time
import zipfile
from CodeFinder import CodeFinder
from MoleculeProcessor import MoleculeProcessor
from NamesProcessor import NamesProcessor
from NMRProcessor import NMRProcessor

@atexit.register
def totalTime():
    print(f"Total execution time: {int(time.time()-startTime)}s")

startTime = time.time()
print("""
BjoStorageSolution
------------------
DEMO
This program will go through the data/ directory, identify and produce a cleaned
version of the data into the <CODE>.zip file
For now the data has to be organized like that:
-/data
 - structure.mol : structure file
 - NAMES.txt : text file containing the names of the molecule/extract
 - CODE.txt : text file containing the code for that molecule
 - <anything>/
      - NMR data

For now it does:
- copy any existing data
- convert the structure to InChI, SMILES and depiction from Mol
- clean up NMR data
- generate a report of the data and structures

""")

class BSSFile:
    code = None
    zipFile = None
    fileName = None
    ontoStore = []

    def setCode(self, code):
        self.code = code

    def __enter__(self):
        fileName = f"./{self.code}.zip"

        self.zipFile = zipfile.ZipFile(fileName, 'w')
        self.fileName = fileName

    def __exit__(self, *args):
        self.zipFile.close()

    def addEntry(self, path, content):
        self.zipFile.writestr(path, content)

    def addFile(self, sourcePath, destinationPath):
        self.zipFile.write(sourcePath, destinationPath)

    def addOnto(self, universal, value):
        self.ontoStore += [{universal: value}]



def nameProcessor(bssFile):
    if not os.path.exists("./data/structure.mol"):
        return False
    return True

# pass1Processors work on the directory and cannot have access to the zipfile yet
pass1Processors = [CodeFinder]
# pass2Processors act inside the zip file, each processor is responsible for adding its own files
pass2Processors = [MoleculeProcessor, NamesProcessor, NMRProcessor]

bssFile = BSSFile()


# TODO: combine these two

for processor in pass1Processors:
    instance = processor(bssFile, "./data")
    ret = instance.process()
    if not ret:
        print(f"Sorry the processor {processor.name} failed.")
        sys.exit(2)


with bssFile:
    for processor in pass2Processors:
        instance = processor(bssFile, "./data")
        ret = instance.process()
        if not ret:
            print(f"Sorry the processor {processor.name} failed.")
            sys.exit(2)
