# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 00:50:19 2016

@author: m.ashenafi
"""

import ImarisLib
import time
import Tkinter, tkFileDialog
import pandas as pd
import math

def DistanceCaluculate(PollenPosition, NucleiPosition):
	DistanceList=[]
	SelectedNuclei=NucleiPosition
	SelectedPollen=PollenPosition
	for b in range(len(SelectedPollen)):
		DistanceX    =   [SelectedPollen[b][0]-c[0] for c in SelectedNuclei]
		DistanceY    =   [SelectedPollen[b][1]-c[1] for c in SelectedNuclei]
		DistanceZ    =   [SelectedPollen[b][2]-c[2] for c in SelectedNuclei]
		Distance     =   [math.sqrt(d**2+e**2+f**2) for d,e,f in zip(DistanceX, DistanceY, DistanceZ)]
		#Distance = [x for x in Distance  if x<10]   
		DistanceList.extend(Distance)
	print DistanceList

def Segment_Surface(vImaris ,vImage, ch, vSFW, name,vLCF, AreaThreshold):
    vROI = None
    vATA = 1
    vATM = 0
    vSFS = '"Area" above '+str(AreaThreshold)+' um^2' if AreaThreshold!=None else AreaThreshold
    vDCI = ch
    vSurface2 = vImaris.GetImageProcessing().DetectSurfaces(vImage, vROI, vDCI,vSFW, vLCF, vATA, vATM, vSFS)
    vSurface2.SetName(name)
    return vSurface2
    
def doIt(id):
	DAPIChannel=2
	vLib = ImarisLib.ImarisLib()
	vImaris = vLib.GetApplication(id)
	vImage=vImaris.GetDataSet()
	vScene                      =    vImaris.GetSurpassScene()
	vFactory		              =	    vImaris.GetFactory()
	GroupOfObjects	           = vFactory.CreateDataContainer()
	GroupOfObjects.SetName('Segmented objects')
	SurfPollen=Segment_Surface(vImaris, vImage, DAPIChannel, 0.8, "Pollen Grain",1, None) # Segment pollen
	#vSigma                      =    0.15
	#vImaris.GetImageProcessing().GaussFilterChannel(vImage,DAPIChannel,vSigma) #Apply Gaussian smoothing filter
	#vSigma=2.4
	#vImaris.GetImageProcessing().SubtractBackgroundChannel(vImage,0,vSigma)
	#aBaseline=float(vImaris.GetDataSet().GetChannelRangeMax(0))/5.0
	#vImaris.GetImageProcessing().BaselineSubtractChannel(vImage,0,aBaseline)
	#SurfNucleus=Segment_Surface(vImaris ,vImage, DAPIChannel, 0.2, "Nuclei",1.15, 1)# Segment nuclei
	#GroupOfObjects.AddChild(SurfNucleus, -1)
	GroupOfObjects.AddChild(SurfPollen, -1)
	vScene.AddChild(GroupOfObjects, -1)
	#NumberOfSurfacesNuclei        =    SurfNucleus.GetNumberOfSurfaces()
	#NumberOfSurfacesPollen        =    SurfPollen.GetNumberOfSurfaces()
	#NucleiPosition=   [SurfNucleus.GetCenterOfMass(surfId)[0] for surfId in range(NumberOfSurfacesNuclei)]
	#PollenPosition =   [SurfPollen.GetCenterOfMass(surfId)[0] for surfId in range(NumberOfSurfacesPollen)]
	#DistanceCaluculate(PollenPosition, NucleiPosition)
	vChannelIndex=2
	vThreshold=-1
	vInside=True
	vImaris.GetImageProcessing().DistanceTransformChannel(vImage,vChannelIndex,vThreshold,vInside)




