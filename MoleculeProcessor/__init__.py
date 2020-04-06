import os
from openbabel import pybel
import tempfile # We need that because openbabel is buggy and doesn't save a png as bytes
from Processor import Processor, ProcessorReturn


class MoleculeProcessor(Processor):
    """Find a molecule.mol file and generate the other parts"""
    name = "MoleculeProcessor"
    help = "Convert the given molecule file to InChI, SMILES and depiction"

    def _process(self):
        molFile = f"{self.directory}/molecule.mol"
        if os.path.exists(molFile):
            print(f"[MoleculeProcessor] Found a molecule. Converting {molFile}")
            mol = list(pybel.readfile("mol", molFile))[0]

            inchikey = mol.write("inchi")
            self.bssFile.addEntry("molecule_inchikey.txt", inchikey)
            self.bssFile.addOnto("http://purl.obolibrary.org/obo/chebi/inchikey", inchikey)

            inchi = mol.write("inchi")
            self.bssFile.addEntry("molecule_inchi.txt", inchi)
            self.bssFile.addOnto("http://purl.obolibrary.org/obo/chebi/inchi", inchi)

            smiles = mol.write("smiles")
            self.bssFile.addEntry("molecule_smiles.txt", smiles)
            self.bssFile.addOnto("http://purl.obolibrary.org/obo/chebi/smiles", smiles)

            svg = mol.write("svg")
            self.bssFile.addEntry("molecule.svg", svg)
            tmp = tempfile.mkstemp(suffix=".png")[1]
            png = mol.draw(False, tmp, True)
            self.bssFile.addFile(tmp, "molecule.png")
            os.remove(tmp)

            self.bssFile.addFile(molFile, "molecule.mol")
            return ProcessorReturn.SUCCESSFUL
        return ProcessorReturn.SKIPPED
