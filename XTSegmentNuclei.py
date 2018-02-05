# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
#   Get surfaces' positions, volumes and intensities
#   Exports result into .csv tables
# Note: This script is calibrated for 3D images DAPI staining of plant nuclei obtained using Leica TCS SP8
# Parameters in XTSegmentNuclei_Parameters.csv file may need to be adjusted for other experiments
# Creator: Mariamawit S. Ashenafi, UZH
# Published on 23.01.2017
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTSegmentNuclei" icon="Python" tooltip="XTSegmentNuclei">
#         <Command>PythonXT::XTSegmentNuclei(%i)</Command>
#       </Item>
#      </Menu>
#    </CustomTools>

import numpy as np
import pandas as pd
import logging
#import tkinter as tk
import tkMessageBox
import os
from Tkinter import *
#from tkinter import filedialog
import tkFileDialog
import time
import datetime
import ImarisLib
#import sys
#sys.path.insert(0, "H:\Python") #I need to add this to import futures on the VMs
#import concurrent.futures

#==============================================================================
# Start of extension
#==============================================================================

#==============================================================================
# Functions required to log and display progress of the plugin
#==============================================================================
def logtime(aTitle):
  # type: (object) -> object
  global gLasttime
  curtime = datetime.datetime.now()
  if (gLasttime is not None):
    diff = (curtime-gLasttime).total_seconds()
  else:
    diff = '-'
  gLasttime = curtime
  print (curtime.ctime(), '[', str(diff), ']', aTitle)
  
#==============================================================================
# Functions required to create new channel to improve segmentation
#==============================================================================

#This function returns the maximum value in a dataframe
def getMax(ListArgument):
    maxList         =   [max(i) for j in ListArgument for i in j]
    maxValue        =   max(maxList)
    return maxValue

#This function returns the image in the form of float value, and get surface table (inside table value=1, outside table value=0)
def Get_Mask_data(Surf):
    vImage			=	vImaris.GetDataSet()
    vImageSizeX 	=	vImage.GetSizeX()
    vImageSizeY 	=	vImage.GetSizeY()
    vImageSizeZ 	=	vImage.GetSizeZ()
    vExtentMinX		=	vImage.GetExtendMinX()
    vExtentMinY		=	vImage.GetExtendMinY()
    vExtentMinZ		=	vImage.GetExtendMinZ()
    vExtentMaxX		=	vImage.GetExtendMaxX()
    vExtentMaxY		=	vImage.GetExtendMaxY()
    vExtentMaxZ		=	vImage.GetExtendMaxZ()
    mask_min		=	[vExtentMinX, vExtentMinY, vExtentMinZ]
    mask_max		=	[vExtentMaxX, vExtentMaxY, vExtentMaxZ]
    mask_size		=	[vImageSizeX, vImageSizeY, vImageSizeZ]
    mask_time		=	0
    mask			=	Surf.GetMask(mask_min[0], mask_min[1], mask_min[2], mask_max[0], mask_max[1], mask_max[2],mask_size[0],mask_size[1], mask_size[2], mask_time)
    mask_values		=	mask.GetDataVolumeFloats(0,0)
    return mask_values

#This function, inside a given surface, selectes voxels containing 50% of the highest intensity, sets there value to the maximum intensity
def SelectVoxels(mask_values, vImage_data, Type, NucSurface, DAPIChannel):
    maxValue=getMax(vImage_data)
#    maxValue= GetSurfIntensity(NucSurface)[DAPIChannel] #get mean value of DAPI inside surface
    Result                  =   []
    for z in range(len(vImage_data)):
        maskY               =   []
        for y in range(len(vImage_data[0])):
            mask_surface	=	mask_values[z][y]
            mask            =	vImage_data[z][y]
#            if Type=="HighIntensity":
#                mask            =   [round(item/maxValue, 0) for item in mask]
#                mask            =	[maxValue if item==1 else item for item in mask]
            mask            =	[x*y for x,y in zip(mask, mask_surface)]
            maskY.append(mask)
        Result.append(maskY)
    return Result

#This function adds a channel in the image, and set the data for this channel
def AddChannel(table, vImage, numberIndex):
    vImage.SetSizeC(numberIndex)
    vImage.SetDataVolumeFloats(table, numberIndex-1, 0)

#==============================================================================
# Functions required to segment DAPI and Immuno channels
#==============================================================================
# Segments the immunostaining channnels:
def GetMaxIntensityOfChannels(vSurpassScene):
    vSurpassScene 		= vImaris.GetSurpassScene()
    vFactory		=	vImaris.GetFactory()
    number_scene_instance	=	range(vSurpassScene.GetNumberOfChildren())
    volume_notFound=True
    i=0
    while i <= number_scene_instance and volume_notFound:
        selection 	=	vSurpassScene.GetChild(i)
        vVolume 	=	vFactory.ToVolume(selection)
        if not vVolume is None:
            volume_notFound=False
        i=i+1
    vAllStatistics 		= 	vVolume.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    Max_intensity_index	= 	[a for a, x in enumerate(vNames) if x == "Data Intensity Max"]
    Max_intensity		=	[vValues[x] for x in Max_intensity_index]
    return Max_intensity

#Segments a channel into surface and returns the surface
def Segment_Surface(vImage, ch, vSFW, name,vLCF, vSFS, vATA, vATM):
    vROI = None
    vDCI = ch
    vSurface2 = vImaris.GetImageProcessing().DetectSurfaces(vImage, vROI, vDCI,vSFW, vLCF, vATA, vATM, vSFS)
    vSurface2.SetName(name)
    return vSurface2

#This function, inside a given surface, selectes voxels containing 50% of the lowest intensity, sets there value to the maximum intensity, and sets the value of the other voxels to 0
def getLowIntensity(mask_values, vImage_data, mask_values2):
    maxValue                =   getMax(vImage_data)
    aVolume                 =   vImaris.GetSurpassSelection()
#    maxValue            =   GetStat(aVolume, "Data Intensity Mean")[DAPIChannel]
#    maxValue            =   GetStat(aVolume, "Data Intensity Max")[DAPIChannel]
    Result                  =   []
    for z in range(len(vImage_data)):
        mask_y              =   []
        for y in range(len(vImage_data[0])):
            mask_surface	=	mask_values[z][y]
            mask_surface2	=	mask_values2[z][y]
            mask            =	vImage_data[z][y]
            # mask=[1 for item==0 else item for item in mask]
            mask            =	[round(item/maxValue, 1) for item in mask]
            mask            =   [0 if item>=0.5 else item for item in mask]
            mask_surface    =   [maxValue if item==0 else 0 for item in mask_surface]
            mask            =	[i*j*k for i,j,k in zip(mask, mask_surface, mask_surface2)]
            mask_y.append(mask)
        Result.append(mask_y)
    return Result

def GetChannelQuantile(vImage_data, Percent):
    IntensityList                  =   []
    for z in range(len(vImage_data)):
        for y in range(len(vImage_data[0])):
            mask            =	vImage_data[z][y]
            IntensityList.extend(mask)
    SegmententationThreshold    =   np.percentile(IntensityList,Percent)
    return SegmententationThreshold

#For a surface containing several surface, this function selectes the most dense surface and creates a surface containing only this selected surface
def SelectSurface(vscene, surf, surfType, groupContainer):
    vTimeIndex          =   0
    vol                 =   GetStat(surf, "Volume")
    if surfType         ==  "Nucleolus":
        area              =   GetStat(surf, "Area")
        volumeToAreaRatio         =   [x/y for x,y in zip(vol, area)]
        SelectedID      =   max(xrange(len(volumeToAreaRatio)), key=volumeToAreaRatio.__getitem__)
    else:
        SelectedID      =   max(xrange(len(vol)), key=vol.__getitem__)
    verts=surf.GetVertices(SelectedID)
    vNormals			=	surf.GetNormals(SelectedID)
    faces=surf.GetTriangles(SelectedID)
    vNucleiSurface		=	vImaris.GetFactory().CreateSurfaces()
    vNucleiSurface.AddSurface(verts, faces, vNormals, vTimeIndex)
    vNucleiSurface.SetName(surfType)
    groupContainer.AddChild(vNucleiSurface, -1)
    vscene.AddChild(groupContainer, -1)
    return vNucleiSurface

#This function first creates a channel by selecting voxel contianing 50% of highest intensity values
def SmoothChannel(chIndex, surf, numberIndex, ParametersList, ChannelName, DAPIChannelColor):
    vImage                  =   vImaris.GetDataSet()
    vScene                  =   vImaris.GetSurpassScene()
    vImage_data             =   vImage.GetDataVolumeFloats(chIndex,0)
    MaskSurface             =   Get_Mask_data(surf)
    result                  =   SelectVoxels(MaskSurface, vImage_data, ChannelName, surf, chIndex)
    AddChannel(result, vImage, numberIndex) #CC channel
    vImage.SetChannelName(numberIndex-1,ChannelName)
    vImage.SetChannelColorRGBA (numberIndex-1, DAPIChannelColor)
    vImaris.GetImageProcessing().GaussFilterChannel(vImage,numberIndex-1,ParametersList[5])
    vImaris.GetImageProcessing().SubtractBackgroundChannel(vImage,numberIndex-1,ParametersList[6])
    
def SegHighIntensity(chIndex, surf, ChannelName, ObjectName, numberIndex, GroupContainer, ParametersList):
    vImage                  =   vImaris.GetDataSet()
#    SmoothChannel(chIndex, surf, numberIndex, ParametersList, ChannelName)
    aSurpassScene           =   vImaris.GetSurpassScene()
#    aVolume                 =   vImaris.GetSurpassSelection()
    vFactory                =   vImaris.GetFactory() 
    numberSceneInstance	    =	  aSurpassScene.GetNumberOfChildren()
    VolumeNotFound          =   True
    i                       =   0
    while i <= numberSceneInstance and VolumeNotFound:
        selection 	          =    aSurpassScene.GetChild(i)
        aVolume                =    vFactory.ToVolume(selection)
        if aVolume is not None:
            VolumeNotFound=False
        i+=1
#    IntensityMax            =   GetStat(aVolume, "Data Intensity Max")[chIndex]
#    IntensityMean            =   GetStat(aVolume, "Data Intensity Mean")[chIndex]
#    IntensitySTD            =   GetStat(aVolume, "Data Intensity StdDev")[DAPIChannel]
#    vSFS                    =   '"Area" above 1.00 um^2'
    vSFS                    =   ''
    vATA                    =   0
#    vATM                    =   IntensityMax/1.7
#    vATM                    =   IntensityMax-(IntensityMax*0.3)
    vImage_data             =   vImage.GetDataVolumeFloats(chIndex,0)
    vATM                    =   GetChannelQuantile(vImage_data, 99)
    result              =   Segment_Surface(vImage, chIndex, 0.0, ObjectName,0.0, vSFS, vATA, vATM)#CC segmented
    return result

#Funtion to get volume and intensity of data item 
def GetStat(surf, FilterString):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    OutputParam	              = 	[float(vValues[a]) for a, x in enumerate(vNames) if x == FilterString]
    return OutputParam

#==============================================================================
# This function: removes all objects created in scene
#==============================================================================
def RemoveObjectsCreated(vScene, ListOfContainers):
    for i in ListOfContainers:
        vScene.RemoveChild(i)

#==============================================================================
# This function:
# Segments chanenls into surface (nucleus, nucleolus and chromocenters) or spots (RNA PolII foci),
# Create masks with the surfaces created
# Count the number of spots in each mask
# And saves results
#==============================================================================
def SegmentAndGetFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, ReturnSurface):
    global vImaris
    global numberIndex
    global DAPIChannel
    vImage              =   	vImaris.GetDataSet()
    if vImage is not None:
        OriginalNumberOfChannel         = vImage.GetSizeC()
        if FileIndex    == 1:
            DAPIChannel= Ask_user(OriginalNumberOfChannel)
        if DAPIChannel<OriginalNumberOfChannel:
            numberIndex     = OriginalNumberOfChannel+1
            ChromocenterSurface =   NucleusSurface  =   NucleolusSurface    =   None
            date			=	str(datetime.datetime.now()).split(" ")[0]
            date			+=	" 00:00:00"
            vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
            vImage                     =    vImaris.GetDataSet()
#==============================================================================
 #        SEGMENT IMAGE
#==============================================================================
            logtime('Segmentation START')
            vScene                     =    vImaris.GetSurpassScene()
            vFactory		           =	vImaris.GetFactory().CreateFactory()
            GroupOfObjects	           =	vFactory.CreateDataContainer()
            GroupOfObjects.SetName('Segmented objects')
            SurfP=Segment_Surface(vImage, DAPIChannel, ParametersList[4], "Rougher Nucleus segmentation", ParametersList[7], '', 1, 0)
            NucleusSurface=SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
            vPathToSaveTables = os.path.join(Result_pathway, "SnapshotNucleus_"+vFileName+".tif")
            vImaris.SaveSnapShot(vPathToSaveTables)
            NucleusSurface.SetVisible(0) #to save snaphot of the the nucleolus and the chromocenters           
            logtime('Nucleus surface segmented - Image '+str(FileIndex))
            if DistanceOptions[1] or DistanceOptions[2]:
                DAPIChannelColor   =   vImage.GetChannelColorRGBA(DAPIChannel)
                SmoothChannel(DAPIChannel, NucleusSurface, numberIndex, ParametersList, "Smoothed DAPI Channel", DAPIChannelColor)
                SmoothedChannel=numberIndex-1
                numberIndex=numberIndex+1
            if DistanceOptions[1]:
                SurfNP=Segment_Surface(vImage, SmoothedChannel, ParametersList[0], "Detailed Nucleus segmentation", ParametersList[7], '', 1, 0)
                if SurfNP.GetNumberOfSurfaces()>0:
                    SurfNP2=Segment_Surface(vImage, SmoothedChannel, ParametersList[1], "Rough Nucleus segmentation",0, '', 1, 0)
                    vImage_data             =   vImage.GetDataVolumeFloats(SmoothedChannel,0)
                    Mask1=Get_Mask_data(SurfNP)
                    Mask2=Get_Mask_data(SurfNP2)
                    resNuc=getLowIntensity(Mask1, vImage_data, Mask2)
                    AddChannel(resNuc, vImage, numberIndex)
                    vImage.SetChannelName(numberIndex-1,"Low DAPI intensity")
                    vImage.SetChannelColorRGBA (numberIndex-1, DAPIChannelColor)
                    NucR=Segment_Surface(vImage, numberIndex-1, ParametersList[2], "Rough Nucleolus segmentation",ParametersList[3], '', 1, 0)
                    numberIndex+=1
                    logtime('Nucleolus surface segmented - Image '+str(FileIndex))
                    NucleolusSurface=SelectSurface(vScene, NucR, "Nucleolus", GroupOfObjects) #Nucleolus is segmented
                    vPathToSaveTables = os.path.join(Result_pathway, "SnapshotNucleolus_"+vFileName+".tif")
                    vImaris.SaveSnapShot(vPathToSaveTables)
                    NucleolusSurface.SetVisible(0) #to save snaphot of the the nucleolus and the chromocenters           

            if DistanceOptions[2]:
                ChromocenterSurface          =  SegHighIntensity(SmoothedChannel, SurfP,"High DAPI intensity",
                                                                 "Chromocenters",
                                                                 numberIndex, GroupOfObjects, ParametersList)#Chromocenters are segmented
                GroupOfObjects.AddChild(ChromocenterSurface, -1)
                vScene.AddChild(GroupOfObjects, -1) 
                vPathToSaveTables = os.path.join(Result_pathway, "SnapshotCC_"+vFileName+".tif")
                vImaris.SaveSnapShot(vPathToSaveTables)
                logtime('Chromocenter surfaces segmented - Image '+str(FileIndex))
            logtime('Segmentation END ')
            if NucleusSurface.GetNumberOfSurfaces()>0:
#==============================================================================
#             GET POSITION OF EVERY SURFACE CREATED
#==============================================================================
                NumberChromocenters     =   ChromocenterSurface.GetNumberOfSurfaces() if DistanceOptions[2] else 0
                SurfLabelA              =   ["Chromocenter"+str(a) for a in range(NumberChromocenters)] if DistanceOptions[2] else []
                SurfLabelB              =   ["Nucleolus"] if DistanceOptions[1] else []
                SurfaceLabel            =   ["Nucleus"] + SurfLabelB + SurfLabelA
                SurfaceLabel            =   pd.DataFrame(SurfaceLabel)
                SurfacePosition         =   [vSurf.GetCenterOfMass(0)[0]  for vSurf in [NucleusSurface, NucleolusSurface] if vSurf is not None]
                SurfaceChromocenter     =   [ChromocenterSurface.GetCenterOfMass(vSurfID)[0] for vSurfID in range(NumberChromocenters)]
                SurfacePosition.extend(SurfaceChromocenter)
                SurfacePosition         =       pd.DataFrame(SurfacePosition)
                logtime('Get surface position for - Image: '+str(FileIndex))
#==============================================================================
#             GET VOLUME OF SURFACES
#==============================================================================
                VolumeNucleolus         =       GetStat(NucleolusSurface, "Volume") if  DistanceOptions[1] else []
                VolumeCC                =       GetStat(ChromocenterSurface, "Volume") if  DistanceOptions[2] else []
                VolumeList              =       GetStat(NucleusSurface, "Volume") + VolumeNucleolus + VolumeCC
                VolumeList              =       pd.DataFrame(VolumeList)
                logtime('Get surface volume')
#==============================================================================
#             GET INTENSITY OF SURFACES
#==============================================================================
                if  DistanceOptions[2]:
                    IntensityListCC         =       GetStat(ChromocenterSurface, "Intensity Mean")
                    IntensityListCC         =       np.array_split(np.array(IntensityListCC),vImage.GetSizeC()) #divide the list into number of channels, one sublist contains intensities for one channel, for surfaces contained in surf
                    IntensityListCC         =       [list(x) for x in IntensityListCC]
                    IntensityListCC         =       pd.DataFrame(IntensityListCC).T
                else: 
                    IntensityListCC         =       pd.DataFrame()
                IntensityListNucleolus      =       GetStat(NucleolusSurface, "Intensity Mean") if  DistanceOptions[1] else []
                IntensitySurface            =       [GetStat(NucleusSurface, "Intensity Mean"),IntensityListNucleolus]
                IntensitySurface            =       pd.concat([pd.DataFrame(IntensitySurface), IntensityListCC])
                IntensitySurface.index      =       range(len(IntensitySurface))
                logtime('Get surface intensity')
#==============================================================================
#             SAVE TABLE
#==============================================================================
                OutPutData              =       pd.concat([SurfaceLabel,SurfacePosition, VolumeList, IntensitySurface], 1)
                IntensityLabel          =       ["IntensityCh"+str(b) for b in range(vImage.GetSizeC())]
                OutPutData.columns      =       ["SurfaceType", "PosX", "PosY", "PosZ", "Volume"]+IntensityLabel
                vPathToSaveTables = os.path.join(Result_pathway, "SurfaceFeatures")
                OutPutData.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.') 
#==============================================================================
#             SAVE MODIFIED IMAGE IF NECESSARY
#==============================================================================        
    #            vPath = os.path.join(Result_pathway, vFileName+".ims")
    #            vImaris.FileSave(vPath, "")
            else:
                print ("There is no channel "+str(DAPIChannel+1)+" for image: "+vFileName)
        else: 
            print ("No nucleus surface detected in image: "+vFileName)
    if BatchProcessing:
        vScene.RemoveChild(GroupOfObjects)
#    os.remove(vFullFileName)
    IsImageCorrect          =     vImage is not None and DAPIChannel<OriginalNumberOfChannel #In both events the images won't be processed
    ReturningVariable=IsImageCorrect if ReturnSurface==False else [IsImageCorrect, NucleusSurface,ChromocenterSurface,  NucleolusSurface, DAPIChannel, GroupOfObjects]
    return ReturningVariable
#==============================================================================
# Functions required to log and display progress of the plugin
#==============================================================================
def logtime(aTitle):
  # type: (object) -> object
  global gLasttime
  curtime = datetime.datetime.now()
  if (gLasttime is not None):
    diff = (curtime-gLasttime).total_seconds()
  else:
    diff = '-'
  gLasttime = curtime
  print (curtime.ctime(), '[', str(diff), ']', aTitle)

#==============================================================================
# Pop-up windows to ask user to set folder pathway and channels to segment
#==============================================================================
class Checkbar(Frame):
  def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
    Frame.__init__(self, parent)
    self.vars = []
    for pick in picks:
      var = IntVar()
      chk = Checkbutton(self, text=pick, variable=var)
      chk.pack(side=side, anchor=anchor, expand=YES)
      self.vars.append(var)
  def state(self):
    return map((lambda var: var.get()), self.vars)

def allstates():
    global User_selection
    global root
    if sum(list(lng.state()))>0:
        User_selection=list(lng.state())
        root.destroy()
    else:
        Message="Please select one of the options."
        Label(root, text=Message).grid(row=1)

def PopUpMessage(OPTIONS, Messge):
  global root
  global lng
  root = Tk()
  label_text	=	Messge
  option		=	OPTIONS
  Label(root, text=label_text).grid(row=0)
  lng = Checkbar(root, option)
  lng.grid(row=2)
  lng.config(relief=GROOVE, bd=2)
  Button(root, text='Quit', fg="red", command=quit).grid(row=4)
  Button(root, text='Submit', fg="darkgreen", command=allstates).grid(row=5)
  root.mainloop()

#This function asks user to set some parameters Immunostaining and DAPI channels.
def Ask_user(numberIndex):
    OPTIONS                         =       range(numberIndex)
    Messge="Please select the DAPI channel: \n Please only choose one channel."
    PopUpMessage(OPTIONS, Messge)
    DAPIChannel	=	[i for i, x in enumerate(User_selection) if x == 1][0]
    return DAPIChannel
  
#==============================================================================
# Function to get parameters for this plugin from file:XTCountSpotPerShell_Parmaeters.csv
#==============================================================================
def GetPluginParameters():
    currentDirectory                         =       os.getcwd()
    AllFilesInDirectory                      =       os.listdir(currentDirectory)
    ParametersList                           =       None
    if "XTCountSpotPerShell_Parameters.csv" in AllFilesInDirectory:
        ParameterData                            =       pd.read_csv("XTSegmentNuclei_Parameters.csv", sep=";", header='infer',decimal='.')
        if "Value" in ParameterData.columns:
            ParametersList                           =       list(ParameterData["Value"])
        else:
            tkMessageBox.showinfo(title="Error", message="Please make sure the 'XTSegmentNuclei_Parameters.csv' file contains a column 'Value' containing the values necessary for this plugin.")
            quit()
    else:
        tkMessageBox.showinfo(title="Error", message="Please make sure there is a 'XTSegmentNuclei_Parameters.csv' in the folder containing the 'XTCountSpotPerShell.py'.")
        quit()
    return ParametersList

#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTSegmentNuclei_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)
        
#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTSegmentNuclei(aImarisId):
	logging.basicConfig(level=logging.DEBUG, filename= "log[XTSegmentNuclei].log")
	try:
		#Declare global variables
		global gLasttime
		global vImaris
		global DAPIChannel
		gLasttime               =    None
		logtime('Extension XTSegmentNuclei START')
		FileNameList            =   []
		#            Step1: Connect to Imaris
		#==============================================================================
		vImarisLib			=	ImarisLib.ImarisLib()
		# Get an imaris object with id aImarisId
		vImaris             =   vImarisLib.GetApplication(aImarisId)
		ParametersList      =   GetPluginParameters()
		# """BEGIN LOOP OVER ALL IMAGES"""
		# Open File and get filename
		if vImaris is not None :
			logtime('Connected to Imaris')
			print("Hello")
			vImaris.GetSurpassCamera().SetPerspective(0) #Set camera to orthographic view 
			ListOfOptions                   =       [["Batch of images", "Just one image"], ["Nucleus", "Nucleolus", "Chromocenters"]]
			ListOfMessages                  =       ["Do you wish to run the script on a batch of images or just on one image already opened?", "Which surface do you wish to create:"]
			UserParameterList               =       []
			for i in range(len(ListOfOptions)):
				OPTIONS                         =       ListOfOptions[i]
				Messge                          =       ListOfMessages[i]
				PopUpMessage(OPTIONS, Messge)
				UserParameterList               =       UserParameterList + [User_selection]
			BatchProcessing	                =	    UserParameterList[0][0]
			DistanceOptions                 =       UserParameterList[1]
			TypeStudy           =   []   # List variable to track the type of study selected by user :  distribution as a function of the nucleus and / or chromocenters and / or nucleolus 
			TypeStudy.append("Nucleus")
			surfaceAnalised     =   [x for x,y in zip(["Chromocenters", "Nucleolus"], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
			SurfaceName         =   '-'.join(surfaceAnalised)
			TypeStudy.extend(SurfaceName)
			ImageIsEmpty        =   False
			FileIndex           =       1 #This  variable is used to count the number of files analysed
			if BatchProcessing  :
				# 		Step2: Here the user is asked to set the path to the folder containing image to be analysed
				#==============================================================================
				root1				    =	        Tk()
				Image_folder			=	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
				root1.destroy()
				Result_pathway          =           os.path.join(Image_folder, "XTSegmentNuclei_Result")
				CreateDirectoryToSaveFiles(Result_pathway)
				AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
				logtime('Get all files')
				AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
				TotalNumberFile         =           len(AllFilesToBeProcessed)
				logtime('Start bach processing')
				if TotalNumberFile  > 0 :
					for vFileName in AllFilesToBeProcessed:
						try:
							vFullFileName = os.path.join(Image_folder, vFileName)
							vImaris.FileOpen(vFullFileName, "")
							logtime('Image processing START image: '+str(vFileName))
							ImageIsEmpty    =           SegmentAndGetFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, False)
							if ImageIsEmpty :
								FileIndex	               +=  		1
								FileNameList.append(vFileName)
							else:
								print ("No image detected in file: "+vFileName)
								logtime('Image processing END')

						except:
								logging.exception("Image="+vFileName+" :")
								print ("There is an issue with file : "+ vFileName)    
								continue #it will not treat this image and move on to the next one                                                      
					df=pd.DataFrame(FileNameList)
					vPathToSaveTables = os.path.join(Result_pathway, "FileName.csv")
					df.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')    
				else:
					tkMessageBox.showinfo(title="Alert", message="There is no .ims or .ics file detected in the selected folder.")
					quit()
			else    :
				TotalNumberFile         =   1
				vFileName               =   vImaris.GetCurrentFileName()
				if vFileName !=""       :
					vFileName           =   vImaris.GetCurrentFileName()
					vFilePath           =   os.path.dirname(vFileName)
					vFullFileName = os.path.join(vFilePath, vFileName)
					Result_pathway      =   os.path.join(vFilePath, "XTSegmentNuclei_Result")
					CreateDirectoryToSaveFiles(Result_pathway)
					vFileName           =   os.path.split(vFileName)[1]
					logtime('Image processing START image: '+str(vFileName))
					ImageIsEmpty        =   SegmentAndGetFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, False)
				if ImageIsEmpty :
					FileIndex	               +=  		1
					print ("NumberOfFile: "+str(FileIndex))
					FileNameList.append(vFileName)					
				else :
					print ("No image is detected. \n Please open an image and select on 'SegmentNuclei' again.")
			print ("All tasks have been completed successfully. \n Resulting files are saved in the folder XTSegmentNuclei_Result")
			raw_input("Press Enter to terminate.")
		else:
			tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
			logtime('Extension XTSegmentNuclei END')
	except:
		logging.exception("Oops:")
        
#==============================================================================
# End of extension
#==============================================================================