#!/usr/bin/env python3
import atexit
import sys
import time

from BSSFile import BSSFile
from BaseProcessors.CodeFinder import CodeFinder
from DataProcessors.MoleculeProcessor import MoleculeProcessor
from DataProcessors.NamesProcessor import NamesProcessor
from DataProcessors.NMRProcessor import NMRProcessor
from FinishingProcessors.ReportGenerator import ReportGenerator


@atexit.register
def total_time():
    print(f"Total execution time: {int(time.time() - startTime)}s")


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

# base processors work on the directory and cannot have access to the zipfile yet
base_processors = [CodeFinder]
# data processors act inside the zip file, each processor is responsible for adding its own files
data_processors = [MoleculeProcessor, NamesProcessor, NMRProcessor]
# finishing processors act inside the zip file, each processor is responsible for adding its own files
finishing_processors = [ReportGenerator]

bssFile = BSSFile()

# TODO: combine these two

for processor in base_processors:
    instance = processor(bssFile, "./data")
    ret = instance.process()
    if not ret:
        print(f"Sorry the processor {processor.name} failed.")
        sys.exit(2)

with bssFile:
    for processor in data_processors:
        instance = processor(bssFile, "./data")
        ret = instance.process()
        if not ret:
            print(f"Sorry the processor {processor.name} failed.")
            sys.exit(2)
    for processor in finishing_processors:
        instance = processor(bssFile, "./data")
        ret = instance.process()
        if not ret:
            print(f"Sorry the processor {processor.name} failed.")
            sys.exit(2)
