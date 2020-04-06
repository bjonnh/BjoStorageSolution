import glob
import os
from openbabel import pybel
from Processor import Processor, ProcessorReturn


class NMRProcessor(Processor):
    """Find NMR datasets and document/clean them"""
    name = "NMRProcessor"
    help = "Process NMR datasets"
    filteredFiles = ["1r", "1i", "2rr", "2ri", "2ir", "2ii"]

    datasets = []

    def findBruker(self):
        return [os.path.dirname(path) for path in glob.glob(f"{self.directory}/**/audita.txt", recursive=True)]

    def processBruker(self, path):
        title = "Untitled (you should give a meaningful title to your NMR experiments)"
        print("  - Adding the experiment (filtering processed data out)")
        for file in glob.glob(f"{path}/**", recursive=True):
            fileName = os.path.basename(file)
            if fileName == "title":
                with open(file, 'r') as f:
                    title = f.read().strip()
            if fileName not in self.filteredFiles:
                newName = self.relPath(file).replace(" ", "_")
                self.bssFile.addFile(file, newName)
        self.datasets += [{"title": title, "path": self.relPath(path).replace(" ", "_"), "type": "bruker"}]

    def _process(self):
        for dataset in self.findBruker():
            print(f"[NMRProcessor] found a Bruker dataset {self.relPath(dataset)}")
            self.processBruker(dataset)
        report = ""
        for dataset in self.datasets:
            report += f"* {dataset['path']} :{dataset['type']}:\n"
            report += "\n".join(["    "+line for line in dataset['title'].split('\n')])
            report += "\n"
        self.bssFile.addEntry("nmr_experiments.txt", report)
        return ProcessorReturn.SKIPPED
