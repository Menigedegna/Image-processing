# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 18:01:55 2017

@author: Nani
"""
%reset
import pandas as pd
import glob, os
import numpy as np


 #    namedf = pd.read_csv(file, skiprows=0, usecols=[1,2,3])

ControlFile=pd.read_csv("Z:/Result0309/AnalyseIntensity/SegmentationCheckObs.csv", skiprows=None)
#OkFiles=ControlFile[FileName[ControlFile$SegmentationNucleus==1]
OkFiles= ControlFile.loc[ControlFile['SegmentationNucleus'] == 1, "FileName"].tolist()
for folderCount in range(2,6): #for the other first concantenate tables in folder then concat all
    os.chdir("Z:/Result0309/s20_Wt/"+str(folderCount)+"/XTSimulateRandomSpots_Result")
    FileName=pd.read_csv("FileName.csv", skiprows=None)["0"].tolist()
    results = pd.DataFrame([])
    for counter, file in enumerate(glob.glob("Intensity_SP*")):
        InputFromFileName=file.split("_")
        FileID=int(InputFromFileName[6].split(".csv")[0])
        #remove files with bad segmentation
        if OkFiles.count(FileName[FileID-1]) == 1: 
            namedf = pd.read_csv(file, skiprows=None)
            pseudo_count = namedf.values[namedf.values > 0].min()/10.
            namedf = np.log10(pseudo_count + namedf)
            namedf = (namedf - namedf.mean(0))/namedf.std(0)
            namedf["SegChannel"]=InputFromFileName[1]
            namedf["Mask"]=InputFromFileName[2]
            namedf["SimulationID"]=InputFromFileName[3]
            namedf["SpotType"]=InputFromFileName[4]
            namedf["NumberSpot"]=InputFromFileName[5]
            results = results.append(namedf, ignore_index=True)
    results.to_csv('Z:/Result0309/s20_Wt/combinedfile_SimIntensity'+str(folderCount)+'.csv')
    
    
results = pd.DataFrame([])
for folderCount in range(1,6): #for the other first concantenate tables in folder then concat all
    namedf = pd.read_csv("Z:/Result0309/s20_Wt/combinedfile_SimIntensity"+folderCount+'.csv', skiprows=None)
    results = results.append(namedf)
results.to_csv('Z:/Result0309/s20_Wt/combinedfile_SimIntensity.csv')
    
    
#    if counter == 3: 
#        break
 

#Check image
from matplotlib import pyplot as plt
from matplotlib import style
import numpy as np
style.use('seaborn-white')

fig,ax = plt.subplots()
ax.set_xscale('log')
ax.set_yscale('log')
ax.scatter(namedf.mean(1),namedf.var(1))
plt.show()

pseudo_count = namedf.values[namedf.values > 0].min()/10.
namedf = np.log10(pseudo_count + namedf)
fig,ax = plt.subplots()
ax.scatter(namedf.mean(1), namedf.std(1))

namedf = (namedf - namedf.mean(0))/namedf.std(0)
fig,ax = plt.subplots()
ax.scatter(namedf.mean(1), namedf.std(1))

from sklearn import metrics
from sklearn import linear_model
from sklearn import cross_validation
from sklearn import ensemble


file="Intensity_SP0_ST1_Sml0_Sim_5000_4.csv"
a=a.split(".csv")[0]
a.split("_")

namedf.head(5)
results.shape
len(results["DAPI"])
#Check data frame type
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
is_string_dtype(results['SpotType'])
is_numeric_dtype(results['DAPI'])

#Apply different functions in columns
for y in agg.columns:
    if (is_string_dtype(agg[y])):
        treat_str(agg[y])
    elif (is_numeric_dtype(agg[y])):
        treat_numeric(agg[y