import os
from openbabel import pybel
from Processor import Processor, ProcessorReturn

class NamesProcessor(Processor):
    """Copy the names"""
    name = "NamesProcessor"
    help = "For now, only copies the NAMES.txt file"

    def _process(self):
        file = f"{self.directory}/NAMES.txt"
        if os.path.exists(file):
            self.bssFile.add_file(file, "NAMES.txt")
            with open(file, 'r') as f:
                self.bssFile.names = [name for name in f.read().split("\n") if name !=""]
            return ProcessorReturn.SUCCESSFUL
        return ProcessorReturn.SKIPPED
