#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 18:04:35 2019

@author: charleshen

@note: Information Content descriptors: http://mordred-descriptor.github.io/documentation/master/api/mordred.InformationContent.html#mordred.InformationContent.InformationContent
"""


from mordred import Calculator, descriptors
import numpy as np

_calc = Calculator(descriptors.InformationContent)
_InfoContentNames = [str(i) for i in _calc.descriptors]


def GetInfoContent(mol):
    """
    #################################################################
    The calculation of InformationContent descriptors (ALL).
    
    Usage:
        
        result=GetInfoContent(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    r = _calc(mol)
    r = r.fill_missing(0)
    return r.asdict()

################################################################
if __name__ =='__main__':

    
    import pandas as pd
    from tqdm import tqdm
    from rdkit import Chem
    
    smis = ['C'*(i+1) for i in range(100)]
    x = []
    for index, smi in tqdm(enumerate(smis), ascii=True):
        m = Chem.MolFromSmiles(smi)
        x.append(GetInfoContent(m))
        
    print(pd.DataFrame(x))  #42