from .fingerprint.atompairs import GetAtomPairFPs
from .fingerprint.avalonfp import GetAvalonFPs
from .fingerprint.rdkitfp import GetRDkitFPs
from .fingerprint.morganfp import GetMorganFPs
from .fingerprint.estatefp import GetEstateFPs
from .fingerprint.maccskeys import GetMACCSFPs
from .fingerprint.pharmErGfp import GetPharmacoErGFPs
from .fingerprint.pharmPointfp import GetPharmacoPFPs
from .fingerprint.pubchemfp import GetPubChemFPs
from .fingerprint.torsions import GetTorsionFPs
from .fingerprint.mhfp6 import GetMHFP6
from .fingerprint.map4 import GetMAP4

from molmap.config import load_config

from rdkit import Chem
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import seaborn as sns

mapfunc = {
           GetMorganFPs:'MorganFP',    
           GetRDkitFPs: 'RDkitFP', 
           GetAtomPairFPs:'AtomPairFP',     
           GetTorsionFPs:'TorsionFP',    
           GetAvalonFPs:'AvalonFP', 
           GetEstateFPs:'EstateFP', 
           GetMACCSFPs:'MACCSFP',  
           GetPharmacoErGFPs:'PharmacoErGFP', 
           GetPharmacoPFPs: 'PharmacoPFP', 
           GetPubChemFPs:'PubChemFP', 
           GetMHFP6:'MHFP6',
           GetMAP4:'MAP4',
          }

mapkey = dict(map(reversed, mapfunc.items()))
colors = sns.palettes.color_palette('hsv', n_colors=len(mapkey)).as_hex()
fps = ['MorganFP','RDkitFP', 'AtomPairFP','TorsionFP',  'AvalonFP','EstateFP','MACCSFP', 'PharmacoErGFP','PharmacoPFP','PubChemFP' ,'MHFP6', 'MAP4']
colormaps = dict(zip(fps, colors))          
colormaps.update({'NaN': '#000000'})

class Extraction:
    
    def __init__(self,  feature_dict = {}):
        """        
        parameters
        -----------------------
        feature_dict: dict parameters for the corresponding fingerprint type, say: {'AtomPairFP':{'nBits':2048}}
        """
        if feature_dict == {}:
            factory = mapkey
            self.flag = 'all'
            cm = colormaps
        else:
            keys = [key for key in set(feature_dict.keys()) & set(mapkey)]
            factory = {}
            cm = {}
            for k, v in mapkey.items():
                if k in keys:
                    factory[k] = mapkey[k]
                    cm[k] = colormaps[k]
            self.flag = 'auto'
        assert factory != {}, 'types of feature %s can be used' % list(mapkey.keys())
            
        self.factory = factory
        self.feature_dict = feature_dict
        _ = self._transform_mol(Chem.MolFromSmiles('CC'))
        self.colormaps = cm        
        self.scaleinfo = load_config('fingerprint', 'scale')
        
    def _transform_mol(self, mol):
        """
        mol: rdkit mol object
        """
        _all = []
        _length = []
        for key,func in self.factory.items():
            kwargs = self.feature_dict.get(key)
            
            if type(kwargs) == dict:
                arr = func(mol, **kwargs)
            else:
                arr = func(mol)
            _length.append(len(arr))
            _all.append(arr)

        concantefp = np.concatenate(_all)
        
        keys = []
        for key, length in zip(self.factory.keys(),  _length):
            keys.extend([(key+str(i), key) for i in range(length)])
            
        bitsinfo = pd.DataFrame(keys, columns=['IDs', 'Subtypes'])
        bitsinfo['colors'] = bitsinfo.Subtypes.map(colormaps)
        self.bitsinfo = bitsinfo            
        return concantefp

    
    def transform(self, smiles):
        '''
        smiles: smile string
        '''
        try:
            mol = Chem.MolFromSmiles(smiles)
            arr = self._transform_mol(mol)
        except:
            #arr = np.nan * np.ones(shape=(len(self.bitsinfo), ))
            arr = np.zeros(shape=(len(self.bitsinfo), ))
            print('error when calculating %s' % smiles)
            
        return arr
    
    
    def batch_transform(self, smiles_list, n_jobs = 1):
        P = Parallel(n_jobs=n_jobs)
        res = P(delayed(self.transform)(smiles) for smiles in tqdm(smiles_list, ascii=True))
        return np.stack(res)
    
def get_fingerprint(input_path):

    all_id = []
    smiles = []

    with open(input_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            drug_id = columns[0]
            smile = columns[1]
            all_id.append(drug_id)
            smiles.append(smile) 

    extraction = Extraction(feature_dict = {})

    smiles_list = smiles

    res = extraction.batch_transform(smiles_list = smiles_list).astype(int)
    drug_fingerprint = pd.DataFrame(res, index=all_id)
    return drug_fingerprint
