# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus into a surface in DAPI channel,
#   Create multiple scaled nucleus surface from nucleus surface to center of mass
#   Get mean and sum DAPI intensity inside each shell (area in-between surfaces)
#   Get volume of shell
#   Exports result into .csv tables
# Note: This script is calibrated for 3D images of plant nuclei conterstained in DAPI obtained using Leica TCS SP8
# Creator: Mariamawit S. Ashenafi, UZH
# Created on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="GetShellIntensity" icon="Python" tooltip="GetShellIntensity">
#         <Command>PythonXT::XTGetShellIntensity(%i)</Command>
#       </Item>
#      </Menu>
#    </CustomTools>

#python library
import numpy as np
import pandas as pd
import logging
import tkMessageBox
import os
from Tkinter import *
import tkFileDialog
import time
import datetime
import itertools
import matplotlib.pyplot as plt
#==============================================================================

#Imaris library
import ImarisLib

#==============================================================================
# Functions required to log and display progress of the plugin
# takes a string as an argument and
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
  print curtime.ctime(), '[', str(diff), ']', aTitle

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

#This function asks user to set some parameters FISH and DAPI channels.
def Ask_user(numberIndex):
    OPTIONS                         =       range(numberIndex)
#    Messge                          =       "Please select the RNA PolII immunostaining channels: \n You can choose multiple channels."
#    PopUpMessage(OPTIONS, Messge)
#    SelectedChanelIndex	            =	    [i for i, x in enumerate(User_selection) if x == 1]
    Messge="Please select the DAPI channel: \n Please only choose one channel."
    PopUpMessage(OPTIONS, Messge)
    DAPIChannel	=	[i for i, x in enumerate(User_selection) if x == 1][0]
#    return SelectedChanelIndex, DAPIChannel
    return  DAPIChannel

#Select voxels with above average intensity -> remove background
def SelectVoxels(surf, mask_values, vImage_data, Type, channel_index):
    meanValue=GetNucleusIntensity(surf, "Mean")[channel_index]
    Result                  =   []
    for z in range(len(vImage_data)):
        maskY               =   []
        for y in range(len(vImage_data[0])):
            mask_surface	=	mask_values[z][y]
            mask            =	vImage_data[z][y]
            mask            =	[x*y for x,y in zip(mask, mask_surface)]
            if Type=="HighIntensity":
                maskNorm            =   [item-meanValue for item in mask]
                mask                =	  [0 if n<=0 else item for n,item in zip(maskNorm, mask)]
            maskY.append(mask)
        Result.append(maskY)
    return Result

#Segments a channel into surface and returns the surface
def Segment_Surface(vImage, ch, vSFW, name,vLCF):
    vROI = None
    vATA = 1
    vATM = 0
    vSFS = ''
    vDCI = ch
    vSurface2 = vImaris.GetImageProcessing().DetectSurfaces(vImage, vROI, vDCI,vSFW, vLCF, vATA, vATM, vSFS)
    vSurface2.SetName(name)
    return vSurface2


#Funtion to get mean intensity inside nucleus surface for all channels selected by the user
def GetNucleusIntensity(surf, TypeOutput):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    IntensitySum        =   [float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity "+ TypeOutput]
    return IntensitySum


#For a surface containing several surface, this function selectes the most dense surface and creates a surface containing only this selected surface
def SelectSurface(vscene, surf, surfType, groupContainer):
    vTimeIndex          =   0
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    volID	            = 	[a for a, x in enumerate(vNames) if x == "Volume"]
    vol		            =	[vValues[x] for x in volID]
    if surfType         ==  "Nucleolus":
        areaID	        = 	[a for a, x in enumerate(vNames) if x == "Area"]
        area		    =	[vValues[x] for x in areaID]
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

#This function returns the volume value of surface
def GetVolume(surf):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    VolumeId	           = 	[a for a, x in enumerate(vNames) if x == "Volume"]
    Volume		           =	float(vValues[VolumeId[0]])
    return Volume


# ==============================================================================
# Functions required to ask user set immuno signal spot diameter for all channels
# ==============================================================================
def fetch(entries, root, ParameterOption):
    global SmothingFactor
    global aLocalContrastFilterWidth
    text = entries.get()
    try:
        x = float(text)
    except ValueError:
        print "You must enter a number in field"
        return
    root.destroy()
    if  ParameterOption=="SmothingSurface":
        SmothingFactor = x
    else :
        aLocalContrastFilterWidth = x


def makeform(root):
    field = "nucleus segmentation"
    row = Frame(root)
    lab = Label(row, width=30, text=field, anchor='w')
    ent = Entry(row)
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent.pack(side=RIGHT, expand=YES, fill=X)
    return ent


def AskUserSmoothingFactor(ParameterOption):
    root = Tk()
    label_text = "Set smooth surface detail " if ParameterOption=="SmothingSurface" else "Set value of Local Contrast Filter Width "
    row = Frame(root)
    lab = Label(row, width=30, text=label_text, anchor='w')
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ents = makeform(root)
    root.bind('<Return>', (lambda event, e=ents: fetch(e, root, ParameterOption)))
    b1 = Button(root, fg="darkgreen", text='Submit',
                command=(lambda e=ents: fetch(e, root, ParameterOption)))
    b1.pack(side=LEFT, padx=5, pady=5)
    b2 = Button(root, fg="darkred", text='Quit', command=root.quit)
    b2.pack(side=LEFT, padx=5, pady=5)
    root.mainloop()



#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTCountSpotPerShell_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)


#==============================================================================
# This function:
# Segments DAPI chanenl into nucleus surface and gets the mean and max DAPI intensity inside nucleus and the volume of the nucleus
# All data is then registered into IntensityMax.csv,IntensityMean.csv and VolumeNuclei.csv
#==============================================================================
def GetImageFeatures(FileIndex, Result_pathway, vFileName, vFullFileName, number_of_shell):
    global vImaris
    global SmothingFactor
    global aLocalContrastFilterWidth
    global IntensityList
    global VolumeList
    global DAPIChannel
    vImage                  =   	vImaris.GetDataSet()
    if vImage is not None:
        numberIndex         =   vImage.GetSizeC()
        if FileIndex        ==  1:
            DAPIChannel     =   Ask_user(numberIndex)
            AskUserSmoothingFactor("SmothingSurface")
            AskUserSmoothingFactor("LocalContrastFilterWidth")
        NucleusSurface      =       None
        date			       =	      str(datetime.datetime.now()).split(" ")[0]
        date			       +=	   " 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        vImage              =       vImaris.GetDataSet()
        logtime('Segmentation START - image_'+str(FileIndex))
        vScene              =       vImaris.GetSurpassScene()
        vFactory		       =	      vImaris.GetFactory().CreateFactory()
        GroupOfObjects	    =    	vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
#        aLocalContrastFilterWidth                =       1.0
        SurfP               =       Segment_Surface(vImage, DAPIChannel, SmothingFactor, "Rougher Nucleus segmentation",aLocalContrastFilterWidth)
        NucleusSurface=SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
        NucleusSurface.SetVisible(False)
        #Normalise image, remove background by setting voxels with below average intensity to 0
        MaskSurface     =   Get_Mask_data(NucleusSurface)
        for channel_index in range(numberIndex):
            vImage_data     =   vImage.GetDataVolumeFloats(channel_index,0)
            CleanedImage    =   SelectVoxels(NucleusSurface, MaskSurface, vImage_data, "HighIntensity", channel_index)
            vImage.SetDataVolumeFloats(CleanedImage, channel_index, 0)#I modify replace the dataset with new dataset
        vPathToSaveTables = os.path.join(Result_pathway, "Snapshot_"+vFileName+".tif")
        vImaris.SaveSnapShot(vPathToSaveTables) #Save snapsht of shells without nucleus
        logtime('Nucleus surface segmented - Image '+str(FileIndex))
        if NucleusSurface is not None:
            numberIndex     =       numberIndex+1
            logtime('Create shell START - Image '+str(FileIndex))
            ListShell=CreateShells(vFactory, NucleusSurface, number_of_shell, GroupOfObjects)
            logtime('Create shell END - Image '+str(FileIndex))
            for shell, nbShell in zip(ListShell,range(number_of_shell, 0, -1)) :
                MaskSurface =   Get_Mask_data(shell)
                AddChannel(MaskSurface, vImage, numberIndex)
                surf        =   Segment_Surface(vImage, numberIndex-1,0.2, "Nucleus shell "+str(nbShell),0) #I create a surface that contains every voxel outside of the chromocenters and inside the nucleus.
#                GroupOfObjects.AddChild(surf, -1)
                Intensity   =   GetNucleusIntensity(surf, "Sum")
                IntensityList.append(Intensity)
                Volume      =   GetVolume(surf)
                VolumeList.append(Volume)
            logtime('Get surface features END - Image '+str(FileIndex))
#        vPath = os.path.join(Result_pathway, vFileName+".ims")
#        vImaris.FileSave(vPath, "") #to save modified file
        vScene.RemoveChild(GroupOfObjects) #remove containter from scene
#        os.remove(vFullFileName) # to delete original file
    return vImage is None

def AddChannel(table, vImage, numberIndex):
    vImage.SetSizeC(numberIndex)
    vImage.SetDataVolumeFloats(table, numberIndex-1, 0)

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

def CreateShells(vFactory, surf, number_of_shell, GroupOfObjects):
    vSurpassScene 		=   vImaris.GetSurpassScene()
    ListOfShells        =   []
    vTimeIndex			=	[0]
    vNormals			=	surf.GetNormals(0)
    vTriangles			=	surf.GetTriangles(0)
    vVertices			=	surf.GetVertices(0)
    vCenterOfMass		=	surf.GetCenterOfMass(0)[0]
    number_of_shell     =   int(number_of_shell)
    NumberOfVerticesPerSurf     =    [len(vVertices)]
    NumberOfTrianglesPerSurf     =    [len(vTriangles)]
    for shell in range(number_of_shell, 0, -1): # for each shell to be created
        NewShellSurface		        =	 vImaris.GetFactory().CreateSurfaces()
        vNewVertices        =   [[(x[0]-vCenterOfMass[0])*(float(shell)/number_of_shell)+vCenterOfMass[0], (x[1]-vCenterOfMass[1])*(float(shell)/number_of_shell)+vCenterOfMass[1], (x[2]-vCenterOfMass[2])*(float(shell)/number_of_shell)+vCenterOfMass[2]] for x in vVertices]
        NewShellSurface.AddSurfacesList(vNewVertices, NumberOfVerticesPerSurf,  vTriangles, NumberOfTrianglesPerSurf,  vNormals, vTimeIndex)
        NewShellSurface.SetName("Shell"+str(shell))
        ListOfShells.append(NewShellSurface)
        GroupOfObjects.AddChild(NewShellSurface, -1)
    return ListOfShells

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTGetShellIntensity(aImarisId):
    logging.basicConfig(level=logging.DEBUG, filename= "log[XTGetShellIntensity].log")
    #Declare global variables
    global gLasttime
    global vImaris
    global IntensityList
    global VolumeList
    gLasttime               =    None
    print "Hello"
    logtime('Extension XTGetShellIntensity START')
    print 'Hello'
    FileNameList            =   []
    VolumeList              =   []
    IntensityList           =   []
    number_of_shell         =   10
    Image_folder = None
    try: #if issue with Imaris connection
#==============================================================================
#                     STEP1 : CONNECT TO IMARIS
#==============================================================================
        vImarisLib			   =	ImarisLib.ImarisLib()
        # Get an imaris object with id aImarisId
        vImaris               =   vImarisLib.GetApplication(aImarisId)
        # """BEGIN LOOP OVER ALL IMAGES"""
        # Open File and get filename
    except:
        logging.exception("Imaris connection issue:")
        print "There was a problem to connect the to Imaris. Please check log file."
        pass
    if vImaris is not None : #Is Imaris connected?
        logtime('Connected to Imaris')
        vImaris.GetSurpassCamera().SetPerspective(0) #Set camera to orthographic view
        FileIndex              =       1 #This  variable is used to count the number of files analysed
#==============================================================================
#              		STEP 2: ASK USER TO SET INPUT PARAMETERS
#==============================================================================
        try: #if issue opening the window to ask user to set parameters
            root1				       =	        Tk()
            Image_folder	       =	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
            root1.destroy()
        except:
            logging.exception("tkFileDialog window issue:")
            print "There was a problem to open window to ask user to set input parameters. Please check log file."
            pass
        if Image_folder is not None:
            Result_pathway          =           os.path.join(Image_folder, "XTGetIntensity_Result")
            CreateDirectoryToSaveFiles(Result_pathway)
            AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
            logtime('Get all files')
            AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
            TotalNumberFile         =           len(AllFilesToBeProcessed)
            logtime('Start bach processing')
            if TotalNumberFile  >   0: #Are there any .ims files?
    #==============================================================================
    #                   STEP 3: OPEN IMAGES IN FOLDER AND START PROCESSING
    #==============================================================================
                for vFileName in AllFilesToBeProcessed:
                    try: #if there is issue with one image (no image inside file or file format issue)
                        vFullFileName   =           os.path.join(Image_folder, vFileName)
                        vImaris.FileOpen(vFullFileName, "")
                        ImageIsEmpty    =           GetImageFeatures(FileIndex, Result_pathway, vFileName, vFullFileName, number_of_shell)
                    except:
                        logging.exception("Image = "+vFileName+" :")
                        print "there is an issue with file : "+ vFileName
                        continue #if there is a problem it moves on to the next file
                    if not ImageIsEmpty :
                        FileIndex	               +=  		1
                        FileNameList.append(vFileName)
    #==============================================================================
    #                     STEP 4: SAVE OUPUT DATA
    #==============================================================================
                vPathToSaveTables       =   os.path.join(Result_pathway, "Result.csv") # pathway to export data
                try:  #if there is issue to format data frame to export output parameters
                    ShellID             =   [range(number_of_shell)]*len(FileNameList)
                    ShellID             =   [x for i in ShellID for x in i]
                    FileNameList        =   [[x]*number_of_shell for x in FileNameList]
                    FileNameList        =   [x for i in FileNameList for x in i]
                    NumberOfChannels    =   max([len(x) for x in IntensityList]) if IntensityList!=[] else [] #get maximum number of channels from all images
                    ResultDf=pd.concat([pd.DataFrame(FileNameList),pd.DataFrame(ShellID), pd.DataFrame(VolumeList), pd.DataFrame(IntensityList)], axis=1)
                    print "Here is the data frame that will be exported:"
                    print ResultDf.loc[0:1,]
                    print "Dimension of the table is : "+str(ResultDf.shape)
                    IntensityLab=["Intensity_Ch"+str(x) for x in range(NumberOfChannels)]
                    ResultDf.columns=["FileName","Shell", "Volume"]+IntensityLab
                    ResultDf.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
                except: #export without formatting
                    logging.exception("Exporting data issue :")
                    print "There was an issue formatting output table: 'Result.csv'. Please check log file."
                    print "The dimension of the table should be columns = 2+number of channels, rows = number of shells * number of images"
                    ResultDf.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
                    pass
            else: #Are there any .ims files?
                tkMessageBox.showinfo(title="Alert", message="There is no .ims or .ics file detected in the selected folder.")
                quit()
        logtime('Extension XTGetShellIntensity END')
        print "Resulting files are saved in the folder XTGetIntensity_Result"
        print "Log file is saved in the same folder containing this plugin"
        raw_input("Press Enter to terminate.")

#==============================================================================
# End of extension
#==============================================================================
