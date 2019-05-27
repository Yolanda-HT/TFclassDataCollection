#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 20:57:10 2019

@author: yolandatiao
"""

########## Match TF family ##########
# Author: Huitian (Yolanda) Diao
# May 26th
# Python 3.7

### Updated on May 24th, 2019
# Loosen motif selection cutoff to include more motifs (Maf, gata...)

########## Import ##########
import os # For changing directory
import csv # For using csv writer
from astropy.io import ascii # For using ascii table to open csv
from astropy.table import Table, Column, join   # For using astropy table functions
import matplotlib.pyplot as plt # For plotting
import pandas as pd # For using pandas in seaborn plot
import numpy as np # For using numpy matrix
import re

########## Self-defined functions ##########
def strInList(strx, listx):
    outlist = []
    for index, i in enumerate(listx):
        if strx in i:
            outlist.append(index)
    return(outlist)
    
def eleInList(elex, listx):
    outlist = []
    outIdx = []
    for i in listx:
        if elex in i:
            outlist.append(i)
            outIdx.append(i.index(elex))
    return([outlist, outIdx])

def lowerList(listx):
    for strx in range(0, len(listx)):
        listx[strx] = listx[strx].lower()
    return(listx)

def ann_ref_TF(inFile):
    inFile = "/Volumes/Yolanda/TFclassDataCollection/MouseTFclassification_noISO.csv"
    outFile = "/Volumes/Yolanda/TFclassDataCollection/MouseTFclassification_noISO_gn_alt.csv"
    refFile = "/Volumes/Yolanda/TFclassDataCollection/Ref/HGNC_symbol_ref_20190527.txt"
    
    refTab = ascii.read(refFile)
    refTab["Synonyms"].fill_value = ""
    refTab["Previous symbols"].fill_value = ""
    ref_symb = list(refTab["Approved symbol"])
    ref_syn = list(refTab["Synonyms"].filled())
    ref_presym = list(refTab["Previous symbols"].filled())
    ref_allsym = [[x]+y.replace(" ","").split(",") + z.replace(" ", "").split(",")  for index, (x, y, z) in enumerate(zip(ref_symb, ref_syn, ref_presym))]
    ref_allsym = [list(filter(None, x)) for x in ref_allsym]
    ref_allsym = [lowerList(x) for x in ref_allsym]
    
    with open(inFile, "r") as fin:
        with open(outFile, "w") as fout:
            rfin = csv.reader(fin, delimiter=",")
            wfout = csv.writer(fout, delimiter=",")
            header = next(rfin)
            wfout.writerow(header + ["Names"])
            for row in rfin:
                row_idx = row[0]
                row_idx_list = row_idx.split(".")
                if len(row_idx_list) == 5:
                    name_x = row[1]
                    name_x_list = name_x.replace("/","").lower().replace(","," ").replace("("," ").replace(")"," ").split(" ")
                    name_x_list = list(filter(None, name_x_list))
                    name_x_sym = []
                    name_x_order = {}
                    name_x_match = {}
                    name_x_count = {}
                    for index, i in enumerate(name_x_list):
                        i_allsym = eleInList(i, ref_allsym)
                        i_allsym_str = [":::".join(sym) for sym in i_allsym[0]]
                        name_x_sym = name_x_sym + i_allsym_str
                        for idx, name in enumerate(i_allsym_str):
                            if name not in name_x_order.keys():
                                name_x_order[name] = index
                                name_x_match[name] = i_allsym[1][idx]
                                name_x_count[name] = 1
                            else:
                                name_x_count[name] = name_x_count[name] + 1
                    name_x_sym = list(set(name_x_sym))
                    
                    # Remove multiple match that are not first order match
                    if len(name_x_sym) > 1:
                        dic_val = name_x_order.values()
                        if min(dic_val) == 0:
                            new_name_x_sym = []
                            for sym in name_x_sym:
                                if name_x_order[sym] == 0:
                                    new_name_x_sym.append(sym)
                            name_x_sym = new_name_x_sym
                    
                    # Remove multiple match that are not first priority name
                    if len(name_x_sym) > 1:
                        dic_val = name_x_match.values()
                        if min(dic_val) == 0:
                            new_name_x_sym = []
                            for sym in name_x_sym:
                                if name_x_match[sym] == 0:
                                    new_name_x_sym.append(sym)
                            name_x_sym = new_name_x_sym
                    
                    # Retry after removing "-"
                    if len(name_x_sym) == 0 or len(name_x_sym) > 1:
                        newname_x = re.sub("^NF-1", "NFI", name_x)
                        newname_x = re.sub("^c-Ets", "Ets", newname_x)
                        newname_x = re.sub("^FIGa", "FIGLA", newname_x)
                        newname_x = re.sub("^L-Myc-1", "Mycl", newname_x)
                        newname_x = re.sub("^C/EBPdelta", "CEBPD", newname_x)
                        newname_x = re.sub("^C/EBPepsilon", "CEBPE", newname_x)
                        newname_x = re.sub("^OLIGO", "OLIG", newname_x)
                        newname_x = re.sub("^IkBdelta", "NFKBID", newname_x)
                        newname_x = re.sub("^IkBepsilon", "NFKBIE", newname_x)
                        newname_x = re.sub("^p63", "tp63", newname_x)
                        newname_x = re.sub("^HOPX-1", "HOPX", newname_x)
                        newname_x = newname_x.replace("-","")
                        name_x_list = newname_x.replace("/","").lower().replace(","," ").replace("("," ").replace(")"," ").split(" ")
                        name_x_list = list(filter(None, name_x_list))
                        name_x_sym = []
                        name_x_order = {}
                        name_x_match = {}
                        name_x_count = {}
                        for index, i in enumerate(name_x_list):
                            i_allsym = eleInList(i, ref_allsym)
                            i_allsym_str = [":::".join(sym) for sym in i_allsym[0]]
                            name_x_sym = name_x_sym + i_allsym_str
                            for idx, name in enumerate(i_allsym_str):
                                if name not in name_x_order.keys():
                                    name_x_order[name] = index
                                    name_x_match[name] = i_allsym[1][idx]
                                    name_x_count[name] = 1
                                else:
                                    name_x_count[name] = name_x_count[name] + 1
                        name_x_sym = list(set(name_x_sym))
                    # Retry after adding "TF"
                    if len(name_x_sym) == 0:
                        newname_x = "TF" + name_x
                        newname_x = newname_x.replace("AP-2delta", "AP2D")
                        newname_x = newname_x.replace("AP-2epsilon", "AP2E")
                        newname_x = newname_x.replace("-", "")
                        name_x_list = newname_x.replace("/","").lower().replace(","," ").replace("("," ").replace(")"," ").split(" ")
                        name_x_list = list(filter(None, name_x_list))
                        name_x_sym = []
                        name_x_order = {}
                        name_x_match = {}
                        name_x_count = {}
                        for index, i in enumerate(name_x_list):
                            i_allsym = eleInList(i, ref_allsym)
                            i_allsym_str = [":::".join(sym) for sym in i_allsym[0]]
                            name_x_sym = name_x_sym + i_allsym_str
                            for idx, name in enumerate(i_allsym_str):
                                if name not in name_x_order.keys():
                                    name_x_order[name] = index
                                    name_x_match[name] = i_allsym[1][idx]
                                    name_x_count[name] = 1
                                else:
                                    name_x_count[name] = name_x_count[name] + 1
                        name_x_sym = list(set(name_x_sym))
                    
                    # Remove multiple match that are not first priority name
                    if len(name_x_sym) > 1:
                        dic_val = name_x_match.values()
                        if min(dic_val) == 0:
                            new_name_x_sym = []
                            for sym in name_x_sym:
                                if name_x_match[sym] == 0:
                                    new_name_x_sym.append(sym)
                            name_x_sym = new_name_x_sym
                    
                    # Remove multiple match that are not first order match
                    if len(name_x_sym) > 1:
                        dic_val = name_x_order.values()
                        if min(dic_val) == 0:
                            new_name_x_sym = []
                            for sym in name_x_sym:
                                if name_x_order[sym] == 0:
                                    new_name_x_sym.append(sym)
                            name_x_sym = new_name_x_sym
                    
                    # Keep the one with largest match count
                    if len(name_x_sym) > 1:
                        dic_val = name_x_count.values()
                        max_count = max(dic_val)
                        new_name_x_sym = []
                        for sym in name_x_sym:
                            if name_x_count[sym] == max_count:
                                new_name_x_sym.append(sym)
                        name_x_sym = new_name_x_sym
                        
                    
                    # Sanity check
                    if len(name_x_sym) == 0:
                        print("No entry found!!! %s"%name_x)   
                        name_x_sym = [name_x]
                    elif len(name_x_sym) > 1:
                        print("Multiple entry found!!! %s"%name_x)
                        
                    
                    newrow = row + [name_x_sym[0].replace(":::", ", ").replace("-", "")]
                    wfout.writerow(newrow)
                else:
                    wfout.writerow(row)
            
    






def anno_family(refFile, inFile):
    inFile = "/Volumes/Yolanda/jycATAC/JYC_DataAnalysis/x_Figures_plottingCodes/0_motifHeatMap/1_combined_motifs/all_sig_motifs.csv"
    outName = inFile.replace(".csv", "_TFfamilyMatched.csv")
    refFile = "/Volumes/Yolanda/jycATAC/JYC_DataAnalysis/z_Info_references/MouseTFclassification_noISO_gn.csv"
    


    

##########---------- Match motifs to families and find consensus motifs
wkdir = "/Volumes/Yolanda/jycATAC/JYC_DataAnalysis/x_Figures_plottingCodes/0_motifHeatMap/1_combined_motifs"
os.chdir(wkdir)







