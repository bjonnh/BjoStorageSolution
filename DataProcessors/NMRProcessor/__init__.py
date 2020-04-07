import glob
import os
from Processor import Processor, ProcessorReturn


class NMRProcessor(Processor):
    """Find NMR datasets and document/clean them"""
    name = "NMRProcessor"
    help = "Process NMR datasets"
    filteredFiles = ["1r", "1i", "2rr", "2ri", "2ir", "2ii"]

    datasets = []

    def find_bruker(self):
        return [os.path.dirname(path) for path in glob.glob(f"{self.directory}/**/audita.txt", recursive=True)]

    def process_bruker(self, path):
        title = "Untitled (you should give a meaningful title to your NMR experiments)"
        print("  - Adding the experiment (filtering processed data out)")
        for file in glob.glob(f"{path}/**", recursive=True):
            file_name = os.path.basename(file)
            if file_name == "title":
                with open(file, 'r') as f:
                    title = f.read().strip()
            if file_name not in self.filteredFiles:
                new_name = self.relPath(file).replace(" ", "_")
                self.bssFile.add_file(file, new_name)
        self.datasets += [{"title": title, "path": self.relPath(path).replace(" ", "_"), "type": "bruker"}]

    def _process(self):
        for dataset in self.find_bruker():
            print(f"[NMRProcessor] found a Bruker dataset {self.relPath(dataset)}")
            self.process_bruker(dataset)
        report = ""
        for dataset in self.datasets:
            self.bssFile.add_dataset_reference(dataset)
            report += f"* {dataset['path']} :{dataset['type']}:\n"
            report += "\n".join(["    "+line for line in dataset['title'].split('\n')])
            report += "\n"
        self.bssFile.add_entry("nmr_experiments.txt", report)
        return ProcessorReturn.SUCCESSFUL
