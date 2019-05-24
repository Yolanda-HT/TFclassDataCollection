#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:29:26 2017

@author: yolandatiao
"""

#####------------------ Import START ------------------#####
import os # For changing directory
import csv # For using csv writer
from astropy.io import ascii # For using ascii table to open csv
from astropy.table import Table, Column    # For using astropy table functions
import glob
import fc_basic_astropy_subprocess as fc # This is a module I have with some useful little functions you can ask me for it
#####------------------ Import END ------------------##### 
 
 
#####------------------ Config START ------------------#####
codedir="/Volumes/Huitian"
wkdir=""
#####------------------ Config END ------------------#####



#####------------------ Self Defined functions START ------------------#####
def tolower(listx):
    outlist=[]
    for i in listx:
        if str(i) != "--":
            outlist.append(i.lower())
        else:
            outlist.append("NaN")
    return outlist
#####------------------ Self Defined functions END ------------------##### 
 
 
#####------------------ Main function START ------------------#####
rnaseqfile="/Volumes/Huitian/Exp174/NormalizedGeneCounts.csv" # RNAseq file name with directory
rtab=ascii.read(rnaseqfile) # Read RNAseq file
#rtab=fc.setcolnames(rtab) # don't use this it's from my module, for solving colunname recognition bug
 
famfile="/Volumes/Huitian/MouseTFclassification/MouseTFclassification_noISO.csv"
famtab=ascii.read(famfile) # Open diction derived transcription family file
famtab=fc.setcolnames(famtab) # don't use this it's from my module, for solving colunname recognition bug
 
rgnlist=list(rtab.columns[0]) # Genenamelist from RNAseq data
famtflist=list(famtab.columns[1]) # Get the column with transcription factors from transcription factor file
famidxlist=list(famtab.columns[0])
 
rgnlistLOW=tolower(rgnlist) # Transform everything to lower case for matcing
famtflistLOW=tolower(famtflist) # Transform everything to lower case for matcing
 
rgnfind=[] # Make an empty list for the results of transcription factor matching
for x in xrange(0,len(famtflistLOW)):
    idx=famidxlist[x]
    if len(idx.split("."))!=5:
        rgnfind.append(None)
    else:
        i=famtflistLOW[x]
        i=i.split("[")[0] # Stripping all the weird characters so that you can match
        i=i.replace(")",",")
        i=i.replace("(",",")
        i=i.replace(" ","")
        i=i.replace("-","")
        i=i.replace("/","")
        ilist=i.split(",") # This is for when there are some variations of names, you will be able to loop match all of them
        ilist=filter(None, ilist)
       
        rgnx=None # Flag for match success
        for n in ilist: # Loop match all variation names
            if n in rgnlistLOW:
                rgnx=n
        rgnfind.append(rgnx) # Add finding results for your match
 
famtab["genename"]=rgnfind # Add match result list to table
ascii.write(famtab,"MouseTFclassification_noISO_GN.csv",format="csv", overwrite=True) # Save

#####------------------ Main function END ------------------#####