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
            self.bssFile.addFile(file, "NAMES.txt")
            return ProcessorReturn.SUCCESSFUL
        return ProcessorReturn.SKIPPED
