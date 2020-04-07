import glob
import re


class CodeFinder:
    """Find a CODE.txt file"""
    name = "CodeFinder"
    help = "Find the CODE.txt file in the data repository"

    def __init__(self, bss_file, directory):
        self.bssFile = bss_file
        self.directory = directory

    def process(self):
        code_file = glob.glob(f"{self.directory}/CODE.txt")
        if len(code_file) == 1:
            with open(code_file[0], 'r') as f:
                code = f.read().strip()
        else:
            return False
        if re.search(r'[^A-Za-z0-9_\-]', code):
            print(f"There are invalid characters in the code")
            return False

        print(f"Found a dataset with the code: {code}")
        self.bssFile.set_code(code)
        return True
