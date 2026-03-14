from .descriptor.autocorr import GetAutocorr, _AutocorrNames
from .descriptor.charge import GetCharge, _ChargeNames
from .descriptor.connectivity import GetConnectivity, _ConnectivityNames
from .descriptor.constitution import GetConstitution, _ConstitutionNames
from .descriptor.estate import GetEstate, _EstateNames
from .descriptor.fragment import GetFragment, _FragmentNames
from .descriptor.kappa import GetKappa, _KappaNames
from .descriptor.moe import GetMOE, _MOENames
from .descriptor.path import GetPath, _PathNames
from .descriptor.property import GetProperty, _PropertyNames
from .descriptor.topology import GetTopology, _TopologyNames
from .descriptor.matrix import GetMatrix, _MatrixNames
from .descriptor.infocontent import GetInfoContent, _InfoContentNames

from molmap.config import load_config

import os
import pandas as pd
import numpy as np
from collections import OrderedDict
from rdkit import Chem
from joblib import Parallel, delayed
from tqdm import tqdm

mapfunc = {GetProperty:'Property', 
           GetConstitution:'Constitution', 
           
           #structure
           GetAutocorr:'Autocorr',
           GetFragment: 'Fragment',
           
           #state
           GetCharge:'Charge',
           GetEstate:'Estate',
           GetMOE:'MOE',
           
           ## graph
           GetConnectivity:'Connectivity', 
           GetTopology:'Topology', 
           GetKappa:'Kappa', 
           GetPath:'Path', 

           GetMatrix:'Matrix', 
           GetInfoContent: 'InfoContent'}

mapkey = dict(map(reversed, mapfunc.items()))

_subclass_ = {'Property':_PropertyNames, 
              'Constitution':_ConstitutionNames,
              'Autocorr':_AutocorrNames,
              'Fragment':_FragmentNames,
              'Charge':_ChargeNames,
              'Estate':_EstateNames,
              'MOE':_MOENames,
              'Connectivity':_ConnectivityNames,
              'Topology':_TopologyNames,
              'Kappa':_KappaNames,
              'Path':_PathNames,
              'Matrix':_MatrixNames,
              'InfoContent':_InfoContentNames}

colormaps = {'Property': '#ff6a00',
             'Constitution': '#ffd500',
             'Autocorr': '#bfff00',
             'Connectivity': '#4fff00',
             'Topology': '#00ff1b',
             'Kappa': '#00ff86', 
             'Path': '#00fff6',             
             'Fragment': '#009eff',
             'Charge': '#0033ff',
             'Estate': '#6568f7',  # #3700ff
             'MOE': '#a700ff',
             'Matrix': '#ff00ed',
             'InfoContent': '#ff0082',
             'NaN': '#000000'}

class Extraction:
    
    def __init__(self,  feature_dict = {}):
        """        
        parameters
        -----------------------
        feature_dict: dict parameters for the corresponding descriptors, say: {'Property':['MolWeight', 'MolSLogP']}
        """
        if feature_dict == {}:
            factory = mapkey
            feature_dict = _subclass_
            self.flag = 'all'
        else:
            factory = {key:mapkey[key] for key in set(feature_dict.keys()) & set(mapkey)}
            feature_dict = feature_dict
            self.flag = 'auto'
        
        assert factory != {}, 'types of feature %s can be used' % list(mapkey.keys())
        self.factory = factory
        self.feature_dict = feature_dict
        keys = []
        for key, lst in self.feature_dict.items():
            if not lst:
                nlst = _subclass_.get(key)
            else:
                nlst = lst
            keys.extend([(v, key) for v in nlst])
        bitsinfo = pd.DataFrame(keys, columns=['IDs', 'Subtypes'])
        bitsinfo['colors'] = bitsinfo.Subtypes.map(colormaps)
        self.bitsinfo = bitsinfo
        self.colormaps = colormaps
        self.scaleinfo = load_config('descriptor','scale')

    def _transform_mol(self, mol):
        """
        mol" rdkit mol object
        """
        _all = OrderedDict()
        
        for key, func in self.factory.items():
            flist = self.feature_dict.get(key)
            dict_res = func(mol)
            if (self.flag == 'all') | (not flist):
                _all.update(dict_res)
            else:
                for k in flist:
                    _all.update({k:dict_res.get(k)})

        arr = np.fromiter(_all.values(), dtype=float)
        arr[np.isinf(arr)] = np.nan #convert inf value with nan            
        return arr
    
    def transform(self, smiles):
        '''
        smiles: smile string
        '''
        try:
            mol = Chem.MolFromSmiles(smiles)
            arr = self._transform_mol(mol)
        except:
            arr = np.nan * np.ones(shape=(len(self.bitsinfo), ))
            print('error when calculating %s' % smiles)
            
        return arr
    
    def batch_transform(self, smiles_list, n_jobs = 1):
        P = Parallel(n_jobs=n_jobs)
        res = P(delayed(self.transform)(smiles) for smiles in tqdm(smiles_list, ascii=True))
        return np.stack(res)
    
def get_descriptor(input_path):

    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)
    txt_file_path = os.path.join(script_dir, "descriptor", "1345_drug_descriptors.txt")
    with open(txt_file_path, 'r', encoding='utf-8') as f:
        descriptors_list = f.read().splitlines()
    descriptors_columns = [int(des) for des in descriptors_list]

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

    res = extraction.batch_transform(smiles_list = smiles_list)
    drug_descriptor = pd.DataFrame(res, index=all_id)[descriptors_columns]
    return drug_descriptor