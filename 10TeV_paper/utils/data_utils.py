import uproot as ur
import awkward as ak
import pathlib

def load_data_from_json(file_path):
    return ak.from_json(pathlib.Path(file_path))


class ROOTDataLoader():
    """
    Note that this currently can only read a single ROOT file at a time.
    """
    def __init__(self):
        self.filename = None
        self.treename = 'ntuple'
        self.file = None
        self.tree = None

    def SetFilename(self,val):
        self.filename = val

    def SetTreename(self,val):
        self.treename=val

    def Load(self):
        if(self.filename is None):
            return
        self.file = ur.open(self.filename)
        self.tree = self.file[self.treename]

    def __getitem__(self,key):
        return self.tree[key].array()
