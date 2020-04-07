import zipfile


class BSSFile:
    __code = None
    __zip_file = None
    __file_name = None
    __dataset_reference_store = []
    __depiction = None
    __ontology_store = []
    __names = []

    def set_code(self, code):
        self.__code = code

    @property
    def code(self):
        """The dataset code"""
        return self.__code

    @property
    def terms(self):
        """The terms associated with this dataset"""
        return self.__ontology_store

    @property
    def datasets(self):
        return self.__dataset_reference_store

    def __enter__(self):
        file_name = f"./{self.__code}.zip"

        self.__zip_file = zipfile.ZipFile(file_name, 'w')
        self.__file_name = file_name

    def __exit__(self, *args):
        self.__zip_file.close()

    def add_dataset_reference(self, dataset_reference):
        self.__dataset_reference_store.append(dataset_reference)

    def add_entry(self, path, content):
        self.__zip_file.writestr(path, content)

    def add_file(self, source_path, destination_path):
        self.__zip_file.write(source_path, destination_path)

    def add_ontological_term(self, universal, value):
        self.__ontology_store += [{universal: value}]

    @property
    def names(self):
        return self.__names

    @names.setter
    def names(self, value):
        self.__names = value

    @property
    def depiction(self):
        return self.__depiction

    @depiction.setter
    def depiction(self, value):
        self.__depiction = value
