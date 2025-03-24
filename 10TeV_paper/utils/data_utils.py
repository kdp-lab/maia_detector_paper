import uproot as ur
import awkward as ak
import pathlib

def load_data_from_json(file_path):
    return ak.from_json(pathlib.Path(file_path))

class JsonDataLoader():
    def __init__(self):
        self.filename = None
        self.dictionary = None
        self.verbose = False

    def SetVerbose(self,val):
        self.verbose = val

    def SetFilename(self,val):
        self.filename = val

    def Load(self):
        if(self.filename is None):
            return
        self.dictionary = ak.from_json(pathlib.Path(self.filename))
        if(self.verbose):
            key = list(self.dictionary.keys())[0]
            nevents = self.dictioanry[key].shape[0]
            print('Loaded {} events.'.format(nevents))

    def __getitem__(self,key):
        return self.dictionary[key]

class ROOTDataLoader():
    # TODO: Having issues when using uproot.open() on a hadd'd file, this causes a crash
    # when opening certain arrays, for example:
    # ValueError: basket 1 in tree/branch /ntuple;1:lc_matched_mcp_pt has the wrong number of bytes (4) for interpretation AsDtype('>f8')
    #
    # It seems that uproot.concatenate() does not have this issue, but I think it has poor memory scaling?
    # NOTE: Do not use uproot.concatenate(), I managed to crash the OSG login node with that...
    #       It really loads all the data into memory, which is exactly what we do not want to do.

    def __init__(self):
        self.filename = None
        self.treename = 'ntuple'
        self.file = None
        self.tree = None
        self.verbose = False

    def SetVerbose(self,val):
        self.verbose = val

    def SetFilename(self,val):
        self.filename = val

    def SetTreename(self,val):
        self.treename=val

    def Load(self):
        if(self.filename is None):
            return
        self.file = ur.open(self.filename)
        self.tree = self.file[self.treename]

        if(self.verbose):
            print('Loaded {} events.'.format(self.tree.num_entries))

    def __getitem__(self,key):
        return self.tree[key].array()

class DataLoader():
    """
    A wrapper class that will open either ROOT ntuples, or JSON files.
    (Please don't use JSON files, they are inefficient and will not scale
    well with memory usage as file sizes get large!)
    """

    def __init__(self):
        self.filename = None
        self.reader = None
        self.mode='ROOT'
        self.verbose=False

    def SetFilename(self,val):
        self.filename = val
        if(self.filename.split('.')[-1].lower() =='json'):
            self.reader = JsonDataLoader()
            self.mode = 'JSON'
        else: # assuming ROOT
            self.reader = ROOTDataLoader()
            self.mode = 'ROOT'
        self.reader.SetFilename(self.filename)
        self.reader.SetVerbose(self.verbose)

    def GetMode(self):
        return self.mode

    def SetVerbose(self,val):
        self.verbose = val

    def Load(self):
        if(self.mode=='JSON'):
            print('Loading JSON data from {}'.format(self.filename))
        elif(self.mode=='ROOT'):
            print('Loading ROOT data from {}:{}'.format(self.filename,self.reader.treename))

        self.reader.Load()

    def __getitem__(self,key):
        return self.reader[key]