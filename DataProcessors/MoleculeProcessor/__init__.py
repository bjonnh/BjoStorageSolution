import os
from openbabel import pybel
import tempfile # We need that because openbabel is buggy and doesn't save a png as bytes
from Processor import Processor, ProcessorReturn


class MoleculeProcessor(Processor):
    """Find a molecule.mol file and generate the other parts"""
    name = "MoleculeProcessor"
    help = "Convert the given molecule file to InChI, SMILES and depiction"

    def _process(self):
        mol_file = f"{self.directory}/molecule.mol"
        if os.path.exists(mol_file):
            print(f"[MoleculeProcessor] Found a molecule. Converting {mol_file}")
            mol = list(pybel.readfile("mol", mol_file))[0]

            inchikey = mol.write("inchikey")
            self.bssFile.add_entry("molecule_inchikey.txt", inchikey)
            self.bssFile.add_ontological_term("http://purl.obolibrary.org/obo/chebi/inchikey", inchikey)

            inchi = mol.write("inchi")
            self.bssFile.add_entry("molecule_inchi.txt", inchi)
            self.bssFile.add_ontological_term("http://purl.obolibrary.org/obo/chebi/inchi", inchi)

            smiles = mol.write("smiles")
            self.bssFile.add_entry("molecule_smiles.txt", smiles)
            self.bssFile.add_ontological_term("http://purl.obolibrary.org/obo/chebi/smiles", smiles)

            svg = mol.write("svg")
            self.bssFile.add_entry("molecule.svg", svg)
            self.bssFile.depiction = "molecule.svg"
            tmp = tempfile.mkstemp(suffix=".png")[1]
            png = mol.draw(False, tmp, True)
            self.bssFile.add_file(tmp, "molecule.png")
            os.remove(tmp)

            self.bssFile.add_file(mol_file, "molecule.mol")
            return ProcessorReturn.SUCCESSFUL
        return ProcessorReturn.SKIPPED
