#!/bin/python3
# Script Name: foreignness.py
# Description: Compute the foreignness (non-selfness) of peptides
# Refer to https://github.com/LukszaLab/NeoantigenEditing

import os
from collections import OrderedDict
from align_neoantigens_to_IEDB import *
from compute_fitness import compute_R

src_dir = os.path.dirname(__file__)

class Foreignness():
    def __init__(self, foreign_file=f'{src_dir}/data/iedb.fasta'):
        # make BLAST db
        if not os.path.isfile(f'{foreign_file}.phr'):
            prepare_blastdb(foreign_file)
        self.db_file = foreign_file
        self.epitopes = load_epitopes(foreign_file)
    
    def __call__(self, peptides, a=22.897590714815188, k=1):
        result_dict = dict()
        # run BLAST
        alignments = run_blastp(peptides, self.db_file, n=100)
        # compute foreignness
        for pept, epitope_ids in alignments.items():
            scores = [align_peptides(pept, self.epitopes[i], blosum62).score for i in epitope_ids] # alignment score
            foreignness = compute_R(scores, a, k) # foreignness
            result_dict[pept] = foreignness
        
        results = list()
        for peptide in peptides:
            results.append(result_dict.get(peptide, 0))
        
        return results