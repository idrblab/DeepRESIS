import sys
import os
import pandas as pd
from Bio import SeqIO
from pathlib import Path
# current_path = Path(__file__).parent.resolve()
# methods_path = current_path / 'methods'
# sys.path.append(str(methods_path))
# prj_path = current_path 
# sys.path.append(prj_path)
# import methods.Methods_all_16_methods as M
from .methods import Methods_all_16_methods as M
current_path = Path(__file__).parent.resolve()

methods = [
    'Open reading frame (1D)',
    'Entropy density of transcript (1D)',
    'Global descriptor (1D)',
    'K-mer (1D)',
    'Codon related (1D)',
    'Pseudo protein related (1D)',
    'Guanine-cytosine related (1D)',
    'Nucleotide related (1D)',
    'EIIP based spectrum (1D)',
    'Solubility lipoaffinity (1D)',
    'Partition coefficient (1D)',
    'Polarizability refractivity (1D)',
    'Hydrogen bond related (1D)',
    'Topological indice (1D)',
    'Molecular fingerprint (1D)'
    ]

dictMe = {
    'Open reading frame (1D)': '1',
    'Entropy density of transcript (1D)': '2_1',
    'Global descriptor (1D)': '2_2',
    'K-mer (1D)': '2_3',
    'Codon related (1D)': '3',
    'Pseudo protein related (1D)': '6',
    'Guanine-cytosine related (1D)': '7',
    'Sequence-intrinsic Features:One-hot encoding (2D)': '13',
    'Sparse encoding (2D)': '15',
    'Structure-based Features:One-hot encoding (2D)':'16',
    'Nucleotide related (1D)':'17',
    'Secondary structure (1D)':'18',
    'EIIP based spectrum (1D)':'18_1',
    'Solubility lipoaffinity (1D)':'19_1',
    'Partition coefficient (1D)':'19_101',
    'Polarizability refractivity (1D)':'19_2',
    'Hydrogen bond related (1D)':'19_3',
    'Topological indice (1D)':'19_4',
    'Molecular fingerprint (1D)':'19_5'
}

def get_ncrna_feature(input_path):

    feature_list = []

    txt_file_path = current_path / 'methods' / 'Data' / 'pending_638_feature.txt'
    with open(txt_file_path, 'r', encoding='utf-8') as f:
        ncrna_feature_name = f.read().splitlines()

    for m in methods:
        T1 = M.switch_meth(dictMe[m], input_path)
        T1.index.name='Seqname'
        feature_list.append(T1)

    combine_feature = pd.concat(feature_list, axis=1)[ncrna_feature_name]
    return combine_feature