from enum import Enum
import os

class ProcessorReturn(Enum):
    SKIPPED = "Skipped"
    FAILED = "Failed"
    SUCCESSFUL = "Successful"

class Processor:
    def __init__(self, bssFile, directory):
        if not hasattr(self, "name"):
            raise ValueError("A Processor needs to have a name defined.")
        self.bssFile = bssFile
        self.directory = directory

    def process(self):
        print(f"[{self.name}] Started")
        ret = self._process()
        print(f"[{self.name}] Finished with status: {ret.value}")
        return ret

    def relPath(self, path):
        """Returns a path without the self.directory"""
        return os.path.relpath(path, start = self.directory)
