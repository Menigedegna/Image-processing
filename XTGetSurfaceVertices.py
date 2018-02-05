# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
#    Export few vertices of surfaces
#   Exports result into .csv tables and save .ims file containing the surfaces, spots created
# Note: This script is calibrated for 3D images of DAPI stained plant nuclei obtained using Leica TCS SP8 or SP5
# Creator: Mariamawit S. Ashenafi, UZH
# Published on 23.01.2018
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTGetSurfaceVertices" icon="Python" tooltip="XTGetSurfaceVertices">
#         <Command>PythonXT::XTGetSurfaceVertices(%i)</Command>
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
def getSurfaceVertices(numberIndex, vSurface, FileIndex, Result_pathway, SurfaceName,GroupOfObjects, NumberOfVertices):
     vFactory		                 =	vImaris.GetFactory().CreateFactory()
     VerticesPosition              =    pd.DataFrame()
     logtime('Get vertices START - image_'+str(FileIndex))
     # 		Quantify spots and intensities
     # ==============================================================================
     for SelectedID in range(vSurface.GetNumberOfSurfaces()):
          vVertices                =    vSurface.GetVertices(SelectedID)
          vVertices                =    pd.DataFrame(vVertices)
          vStep                    =    len(vVertices)/NumberOfVertices
          if vStep>2: 
               SelectIndex              =    range(0, len(vVertices), vStep)
               vVertices                =    vVertices.iloc[SelectIndex,]
          verts                    =    vVertices.values.tolist()
          CreateSpots(vFactory,verts , SurfaceName+"Vertices", GroupOfObjects, 0.025)
          VerticesPosition         =    VerticesPosition.append(vVertices)
     #==============================================================================
     #                 Export spot positions
     #==============================================================================
     VerticesPosition.index        =    range(len(VerticesPosition))
     vPathToSaveTables = os.path.join(Result_pathway,SurfaceName+"Vertices")
     VerticesPosition.to_csv(path_or_buf=vPathToSaveTables+"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
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
    return NucleusSurface, NucleolusSurface, ChromocenterSurface, ImmunoSpotNames, ImmunoSpotList

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
        numberIndex     +=1
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        GroupOfObjects	           =	vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
#==============================================================================
# SEGMENT SURFACES
#==============================================================================
#==============================================================================
#         Reset data coordinate so that the min X,Y,Z=0, this is necessary to simulate random spots
#==============================================================================
        vExtentMinX		=	vImage.GetExtendMinX()
        vExtentMinY		=	vImage.GetExtendMinY()
        vExtentMinZ		=	vImage.GetExtendMinZ()
        vExtentMaxX		=	vImage.GetExtendMaxX()
        vExtentMaxY		=	vImage.GetExtendMaxY()
        vExtentMaxZ		=	vImage.GetExtendMaxZ()
        vImage.SetExtendMaxX(vExtentMaxX-vExtentMinX) 
        vImage.SetExtendMaxY(vExtentMaxY-vExtentMinY) 
        vImage.SetExtendMaxZ(vExtentMaxZ-vExtentMinZ) 
        vImage.SetExtendMinX(0.00) 
        vImage.SetExtendMinY(0.00) 
        vImage.SetExtendMinZ(0.00)
        if DoSegmentation:
            IsImageCorrect, NucleusSurface,ChromocenterSurface,  NucleolusSurface, DAPIChannel, GroupOfObjects=SN.SegmentAndGetFeatures(vImage, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, True)
            logtime('Segmentation END - image_'+str(FileIndex))
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface = GetSegmentedObjects(DistanceOptions, FileIndex)
#==============================================================================
# GET VERTICES
#==============================================================================
        if DistanceOptions[0] and NucleusSurface is not None :
             getSurfaceVertices(numberIndex, NucleusSurface, FileIndex, Result_pathway, "Nucleus ", GroupOfObjects, 50000)
        if DistanceOptions[2] and ChromocenterSurface is not None :
             getSurfaceVertices(numberIndex, ChromocenterSurface, FileIndex, Result_pathway, "Chromocenters ",GroupOfObjects,  500)
        if DistanceOptions[1] and NucleolusSurface is not None :
             getSurfaceVertices(numberIndex, NucleolusSurface, FileIndex, Result_pathway, "Nucleolus ",GroupOfObjects, 10000)
        ListOfContainers.append(GroupOfObjects)
    else:
        print ("No image detected in file: "+vFileName)
        quit()
    if len(ListOfContainers)>0 and BatchProcessing:
        RemoveObjectsCreated(vScene, ListOfContainers)
#    os.remove(vFullFileName)
    return vImage is None

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
    if "XTCountSpotPerShell_Parameters.csv" in AllFilesInDirectory:
        ParameterData                            =       pd.read_csv("XTCountSpotPerShell_Parameters.csv", sep=";", header='infer',decimal='.')
        if "Value" in ParameterData.columns:
            ParametersList                           =       list(ParameterData["Value"])
        else:
            tkMessageBox.showinfo(title="Error", message="Please make sure the 'XTCountSpotPerShell_Parameters.csv' file contains a column 'Value' containing the values necessary for this plugin.")
            quit()
    else:
        tkMessageBox.showinfo(title="Error", message="Please make sure there is a 'XTCountSpotPerShell_Parameters.csv' in the folder containing the 'XTCountSpotPerShell.py'.")
        quit()
    return ParametersList

#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTGetSurfaceVertices_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTGetSurfaceVertices(aImarisId):
	logging.basicConfig(level=logging.DEBUG, filename= "log[XTGetSurfaceVertices].log")
	try:
		#Declare global variables
		global gLasttime
		global vImaris
		global SelectedChanelIndex
		gLasttime               =    None
		logtime('Extension XTGetSurfaceVertices START')
		print ("Hello")
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
			vImaris.GetSurpassCamera().SetOrthographic(True) #Set camera to orthographic view 
			vImaris.GetSurpassCamera().Fit() #Sets the zoom and the position so that the bounding box of all visible objects fits into the window 
			ListOfOptions                   =       [["Batch of images", "Just one image"], ["Segment & Get vertices", "Get vertices"], ["Nucleus", "Nucleolus", "Chromocenters"]]
			ListOfMessages                  =       ["Do you wish to run the script on a batch of images or just on one image already opened?", "Do you wish to do automated segmentation?", "Do you wish to analyse RNA PolII distribution's as a function of the:"]
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
				Result_pathway          =   os.path.join(r"Y:\Result0309\test", FolderName, "XTGetSurfaceVertices_Result")
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
        						ImageIsEmpty    =           GetImageFeatures(FileIndex, Result_pathway, vFileName, DoSegmentation, ParametersList, BatchProcessing, vFullFileName, DistanceOptions)
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
					Result_pathway      =   os.path.join(vFilePath, "XTGetSurfaceVertices_Result")
					CreateDirectoryToSaveFiles(Result_pathway)
					vFileName           =   os.path.split(vFileName)[1]
					ImageIsEmpty        =   GetImageFeatures(FileIndex, Result_pathway, vFileName, DoSegmentation, ParametersList, BatchProcessing, vFullFileName, DistanceOptions)
				if not ImageIsEmpty :
					FileIndex	               +=  		1
					FileNameList.append(vFileName)					
				if ImageIsEmpty :
					tkMessageBox.showinfo(title="Alert", message="No image is detected. \n Please open an image and select on 'CountSpotPerShell' again.")
					quit()
			logtime('XTGetSurfaceVertices extension done')
			print ("All tasks have been completed successfully. \n Resulting files are saved in the folder XTGetSurfaceVertices_Result")
			raw_input("Press Enter to terminate.")
		else:
			tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
			logtime('Extension XTGetSurfaceVertices END')
	except:
		logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
