# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
#    Segment FISH or Immunostaining signal into spots or surfaces
#   Exports items position and features into a .csv tables and save .ims file 
# Note: This script is designed for whole mount immunostaining, FISH of plant tissue conterstained in DAPI. 3D images are acquired using high resolution microscopy. Parameters in XTWholeMountSegmentation_Parameters.csv need to be adjusted to resolution of image.
# Creator: Mariamawit S. Ashenafi, UZH
# Created on 11.04.2018
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTWholeMountSegmentation" icon="Python" tooltip="XTWholeMountSegmentation">
#         <Command>PythonXT::XTWholeMountSegmentation(%i)</Command>
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
import sys
sys.path.insert(0, "H:\Python") #I need to add this to import futures on the VMs
import concurrent.futures
import XTSegmentNuclei as SN

#==============================================================================
# Start of extension
#==============================================================================

# Function to create spots
def CreateSpots(vFactory, aPositionsXYZ, SpotName, groupContainer, aRadius):
    vSpot			        =	vFactory.CreateSpots()
    aIndicesT				=	[0.0]*len(aPositionsXYZ)
    aRadii					=	[aRadius]*len(aPositionsXYZ)
    vSpot.Set(aPositionsXYZ,aIndicesT,aRadii)
    vSpot.SetName(SpotName)
    groupContainer.AddChild(vSpot, -1)
    return vSpot



# Function to quantify spots in image
def getSurfaceVertices(numberIndex, vSurface, FileIndex, Result_pathway, SurfaceName,GroupOfObjects, LayerSize):
     vFactory		                 =	vImaris.GetFactory().CreateFactory()
     NumberOfSurfaces = vSurface.GetNumberOfSurfaces()
     logtime('Get vertices START - image_'+str(FileIndex))
#==============================================================================
#  Get vertices per suface and add a layer for each nucleus
# 
#==============================================================================
#     vVertices=[]
#     vNumberOfVerticesPerSurface=[]
#     vTriangles=[]
#     vNumberOfTrianglesPerSurface=[]
#     vNormals=[]
     NewShellSurface		        =	 vImaris.GetFactory().CreateSurfaces()
     for SelectedID in range(NumberOfSurfaces):
          aVertices                =  vSurface.GetVertices(SelectedID)
#          vNumberOfVerticesPerSurface=vNumberOfVerticesPerSurface+[len(aVertices)]
          vTriangles                =  vSurface.GetTriangles(SelectedID)
#          vTriangles.extend(aTriangles)
#          vNumberOfTrianglesPerSurface=vNumberOfTrianglesPerSurface+[len(aTriangles)]
          vNormals                =  vSurface.GetNormals(SelectedID)
#          vNormals.extend(aNormals)
          vCenterOfMass		=	vSurface.GetCenterOfMass(SelectedID)[0]
#          aVertices                =    pd.DataFrame(aVertices)
#          vStep                    =    len(vVertices)/NumberOfVertices
#          if vStep>2: 
#               SelectIndex              =    range(0, len(vVertices), vStep)
#               vVertices                =    vVertices.iloc[SelectIndex,]
#          LayerMat=[x-y for x,y in zip(LayerMat, vCenterOfMass)]
#          LayerMat=pd.DataFrame(LayerMat).T
#          LayerMat = pd.concat([LayerMat]*len(aVertices))
#          aVertices=aVertices.add(LayerMat)
#          vVertices                    =    aVertices.values.tolist()
          vVertices=[[x[0]-LayerSize if x[0]<vCenterOfMass[0] else x[0]+LayerSize, x[1]-LayerSize if x[1]<vCenterOfMass[1] else x[1]+LayerSize, x[2]-LayerSize if x[2]<vCenterOfMass[2] else x[2]+LayerSize] for x in aVertices]
#          vVertices        =   [[x[0]+LayerSize, x[1]+LayerSize, x[2]+LayerSize] for x in aVertices]
          vTimeIndexPerSurface           =  0
          NewShellSurface.AddSurface (vVertices, vTriangles,vNormals,vTimeIndexPerSurface)
#     NewShellSurface.AddSurfacesList(vVertices,vNumberOfVerticesPerSurface,vTriangles,vNumberOfTrianglesPerSurface,vNormals,vTimeIndexPerSurface)
     NewShellSurface.SetName("Adjusted nuclei")
     GroupOfObjects.AddChild(NewShellSurface, -1)
     logtime('Get vertices END - image_'+str(FileIndex))


#==============================================================================
# This function: removes all objects created in scene
#==============================================================================
def RemoveObjectsCreated(vScene, ListOfContainers):
    for i in ListOfContainers:
        vScene.RemoveChild(i)

def GetSegmentedObjects(DistanceOptions,FileIndex):
    logtime('Object detection START - image_'+str(FileIndex))
    NucleusSurface  =   NucleolusSurface    =   ChromocenterSurface =   None
    vScene                     =    vImaris.GetSurpassScene()
    vFactory		           =    vImaris.GetFactory()
    numberSceneInstance	       =	vScene.GetNumberOfChildren()
    ContainerNotFound          =    True
    i                          =    0
    while i <= numberSceneInstance and ContainerNotFound:
        selection 	           =	vScene.GetChild(i)
        vContainer             =	vFactory.ToDataContainer(selection)
        if vContainer is not None:
            ContainerNotFound      =    False
        i=i+1
    ContainerName             =    vContainer.GetName()
    if vContainer is not None and ContainerName=="Segmented objects":
        numberSceneInstance  =    vContainer.GetNumberOfChildren()
        i                    =    0
#        ImmunoSpotNames        =    []
#        ImmunoSpotList         =    []
        while i <= numberSceneInstance :
            selection 	          =    vContainer.GetChild(i)
            vObject               =    vFactory.ToSurfaces(selection)
            if vObject is not None:
                if vObject.GetName() == "Nucleus" and DistanceOptions[0]        :   NucleusSurface      = vObject
                if vObject.GetName() == "Nucleolus" and DistanceOptions[1]      :   NucleolusSurface    = vObject
                if vObject.GetName() == "Chromocenters" and DistanceOptions[2]  :   ChromocenterSurface = vObject
#            vObject               =    vFactory.ToSpots(selection)
#            if vObject is not None:
#                PositionTable =     vObject.GetPositionsXYZ()
#                if "Immuno" in vObject.GetName() and PositionTable is not None:
#                    ImmunoSpotNames     = ImmunoSpotNames   +    [vObject.GetName()]
#                    ImmunoSpotList      = ImmunoSpotList    +    [vObject]
            i+=1
        logtime('Object detection END - image_'+str(FileIndex))
#    return NucleusSurface, NucleolusSurface, ChromocenterSurface, ImmunoSpotNames, ImmunoSpotList
    return NucleusSurface, NucleolusSurface, ChromocenterSurface

#==============================================================================
# This function:
# Segments chanenls into surface (nucleus, nucleolus and chromocenters) or spots (RNA PolII foci),
# Create masks with the surfaces created
# Count the number of spots in each mask
# And saves results
#==============================================================================
def GetImageFeatures(FileIndex, Result_pathway, vFileName, DoSegmentation, ParametersList, BatchProcessing, vFullFileName, DistanceOptions):
    global vImaris
    global numberIndex
    global DAPIChannel
    global ListOfContainers
    ListOfContainers    =   []   # To keep track of all the containers I will create in scene, so that I can remove them after saving image
    ChromocenterSurface =   NucleusSurface  =   NucleolusSurface    =   None
    vImage              =   	vImaris.GetDataSet()
    vScene                     =    vImaris.GetSurpassScene()
    vFactory		           =	vImaris.GetFactory().CreateFactory()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        #Ask user to set FISH channel to segment into surface
    if FileIndex    == 1:
        FISHChannelList= SN.Ask_user(numberIndex, "FISH")
        numberIndex     +=1
        GroupOfObjects	           =	vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
#==============================================================================
# SEGMENT SURFACES
#==============================================================================
        if DoSegmentation:
            IsImageCorrect, NucleusSurface,ChromocenterSurface,  NucleolusSurface, DAPIChannel, GroupOfObjects=SN.SegmentAndGetFeatures(vImage, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, True)
            logtime('DAPI Segmentation END - image_'+str(FileIndex))
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface = GetSegmentedObjects(DistanceOptions, FileIndex)
#Adjust nucleus surfaces
        if ParametersList[8]>0:
             getSurfaceVertices(numberIndex, NucleusSurface, FileIndex, Result_pathway, "Nucleus ", GroupOfObjects, ParametersList[8])
#Segment FISH channel into surfaces
        for FISHChannel in FISHChannelList:
             fishChannelColor   =   vImage.GetChannelColorRGBA(FISHChannel)
             SN.SmoothChannel(FISHChannel, NucleusSurface, numberIndex, ParametersList, "Smoothed FISH Channel", fishChannelColor)
             SmoothedChannel=numberIndex-1
             FISHSurface          =  SN.SegHighIntensity(SmoothedChannel, NucleusSurface,"High FISH intensity","FISH Surface",numberIndex, GroupOfObjects, ParametersList)
             GroupOfObjects.AddChild(FISHSurface, -1)
             vScene.AddChild(GroupOfObjects, -1) 
             logtime('FISH Channel-'+str(FISHChannel)+' Segmentation END - image_'+str(FileIndex))
             ListOfContainers.append(GroupOfObjects)
             ResultFileName="FISHCh"+str(FISHChannel)+"_SurfaceFeatures"
             SN.ExtractSurfaceFeatures([FISHSurface],["FISH"], Result_pathway, FileIndex, vImage, ResultFileName, NucleusSurface)
             logtime('FISH Channel-'+str(FISHChannel)+' surface features saved END - image_'+str(FileIndex))
        vPath = os.path.join(Result_pathway, vFileName+".ims")
        vImaris.FileSave(vPath, "")
    else:
        print ("No image detected in file: "+vFileName)
        quit()
    if len(ListOfContainers)>0 and BatchProcessing:
        RemoveObjectsCreated(vScene, ListOfContainers)
 #    os.remove(vFullFileName)
    return vImage is None
 
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
  
#==============================================================================
# Function to get parameters for this plugin from file:XTCountSpotPerShell_Parmaeters.csv
#==============================================================================
def GetPluginParameters():
    currentDirectory                         =       os.getcwd()
    AllFilesInDirectory                      =       os.listdir(currentDirectory)
    ParametersList                           =       None
    if "XTWholeMountSegmentation_Parameters.csv" in AllFilesInDirectory:
        ParameterData                            =       pd.read_csv("XTWholeMountSegmentation_Parameters.csv", sep=";", header='infer',decimal='.')
        if "Value" in ParameterData.columns:
            ParametersList                           =       list(ParameterData["Value"])
        else:
            tkMessageBox.showinfo(title="Error", message="Please make sure the 'XTWholeMountSegmentation_Parameters.csv' file contains a column 'Value' containing the values necessary for this plugin.")
            quit()
    else:
        tkMessageBox.showinfo(title="Error", message="Please make sure there is a 'XTWholeMountSegmentation_Parameters.csv' in the folder containing the 'XTCountSpotPerShell.py'.")
        quit()
    return ParametersList

#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTWholeMountSegmentation_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTWholeMountSegmentation(aImarisId):
	logging.basicConfig(level=logging.DEBUG, filename= "log[XTWholeMountSegmentation].log")
	try:
		#Declare global variables
		global gLasttime
		global vImaris
		global SelectedChanelIndex
		gLasttime               =    None
		FISHChannel           =    None
		logtime('Extension XTWholeMountSegmentation START')
		print ("Hello!")
		FileNameList            =   []
		#            Step1: Connect to Imaris
		#==============================================================================
		vImarisLib			=	ImarisLib.ImarisLib()
		# Get an imaris object with id aImarisId
		vImaris             =   vImarisLib.GetApplication(aImarisId)
		ParametersList      =   GetPluginParameters()
		SN.vImaris=vImaris
		SN.gLasttime=gLasttime
		# """BEGIN LOOP OVER ALL IMAGES"""
		# Open File and get filename
		if vImaris is not None :
			logtime('Connected to Imaris')
#Set camera to orthographic view 
			vImaris.GetSurpassCamera().SetOrthographic(True) 
			vImaris.GetSurpassCamera().Fit() #Sets the zoom and the position so that the bounding box of all visible objects fits into the window 
			ListOfOptions                   =       [["Batch of images", "Just one image"], ["Segment & Get Features", "Get Features"], ["Nucleus", "Nucleolus", "Chromocenters"]]
			ListOfMessages                  =       ["Do you wish to run the script on a batch of images or just on one image already opened?", "Do you wish to do automated segmentation?", "Which nucleus features do you wish to segment?"]
			UserParameterList               =       []
			for i in range(len(ListOfOptions)):
				OPTIONS                         =       ListOfOptions[i]
				Messge                          =       ListOfMessages[i]
				PopUpMessage(OPTIONS, Messge)
				UserParameterList               =       UserParameterList + [User_selection]
			BatchProcessing	                =	    UserParameterList[0][0]
			DoSegmentation	                =	    UserParameterList[1][0]
			DistanceOptions                 =       UserParameterList[2]
			FileIndex                       =       1 #This  variable is used to count the number of files analysed
			if BatchProcessing  :
				# 		Step2: Here the user is asked to set the path to the folder containing image to be analysed
				#==============================================================================
				root1				    =	        Tk()
				Image_folder			=	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
				root1.destroy()
				FolderName          =   os.path.basename(Image_folder)
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
        #						with concurrent.futures.ProcessPoolExecutor() as executor:
        						ImageIsEmpty    =           GetImageFeatures(FileIndex, Result_pathway, vFileName, DoSegmentation, ParametersList, BatchProcessing, vFullFileName, DistanceOptions, FISHChannel)
        						if not ImageIsEmpty :
        							FileIndex	               +=  		1
        							FileNameList.append(vFileName)
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
				ImageIsEmpty            =   True
				if vFileName !=""       :
					vFileName           =   vImaris.GetCurrentFileName()
					vFilePath           =   os.path.dirname(vFileName)
					vFullFileName = os.path.join(vFilePath, vFileName)
					Result_pathway      =   os.path.join(vFilePath, "XTWholeMountSegmentation_Result")
					CreateDirectoryToSaveFiles(Result_pathway)
					vFileName           =   os.path.split(vFileName)[1]
					ImageIsEmpty        =   GetImageFeatures(FileIndex, Result_pathway, vFileName, DoSegmentation, ParametersList, BatchProcessing, vFullFileName, DistanceOptions)
				if not ImageIsEmpty :
					FileIndex	               +=  		1
					FileNameList.append(vFileName)					
				if ImageIsEmpty :
					tkMessageBox.showinfo(title="Alert", message="No image is detected. \n Please open an image and select on 'CountSpotPerShell' again.")
					quit()
			logtime('XTWholeMountSegmentation extension done')
			print ("All tasks have been completed successfully. \n Resulting files are saved in the folder XTWholeMountSegmentation_Result")
			raw_input("Press Enter to terminate.")
		else:
			tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
			logtime('Extension XTWholeMountSegmentation END')
	except:
		logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
