#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 20:29:36 2019

@author: charleshen

@note: calculation of 606 molecular autocorrelation  descriptors, with three types: MoreauBroto, Moran, Geary, 

ref 0: Moriwaki, Hirotomo, et al. "Mordred: a molecular descriptor calculator." Journal of cheminformatics 10.1 (2018): 4.
ref 1: http://www.rguha.net/writing/notes/desc/node2.html
ref 2: Todeschini and Consoni “Descriptors from Molecular Geometry” Handbook of Chemoinformatics http://dx.doi.org/10.1002/9783527618279.ch37

"""

from mordred import Calculator, descriptors

import numpy as np

_calc = Calculator(descriptors.Autocorrelation)

_AutocorrNames = [str(i) for i in _calc.descriptors]


def GetAutocorr(mol):
    r = _calc(mol)
    r = r.fill_missing(0)
    return r.asdict()

if __name__ =='__main__':
    
    import pandas as pd
    from tqdm import tqdm
    from rdkit import Chem
    smis = ['C'*(i+1) for i in range(100)]
    x = []
    for index, smi in tqdm(enumerate(smis), ascii=True):
        m = Chem.MolFromSmiles(smi)
        x.append(GetAutocorr(m))
        
    print(pd.DataFrame(x))  #606 (sparse maybe)