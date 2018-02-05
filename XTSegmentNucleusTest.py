# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:08:35 2017

@author: Nani
"""
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
#   Create masks using these surfaces
#   Segments RNA PolII immunostaining signal into spots in Immunostaining channels,
#   Count the number of spots in each shell created
#   Get intensity for selected channles for all segmentted spots
#   Exports result into .csv tables and save .ims file containing the surfaces, spots and masks created
# Note: This script is calibrated for 3D images of RNA PolII immunostaining in plant nuclei obtained using Leica TCS SP8
# Creator: Mariamawit S. Ashenafi, Célia Baroux, UZH
# Published on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="SegmentNucleusTest" icon="Python" tooltip="SegmentNucleusTest">
#         <Command>PythonXT::XTSegmentNucleusTest(%i)</Command>
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

#==============================================================================
# Start of extension
#==============================================================================
#This function returns the volume value of surface
def GetNucParameters(surf):
    vAllStatistics              =     surf.GetStatistics()
    vNames                      =     vAllStatistics.mNames
    vValues                     =     vAllStatistics.mValues
    SelectedId                  =     [a for a, x in enumerate(vNames) if x == "Volume"]
    Volume                      =     [float(vValues[id]) for id in SelectedId] if len(SelectedId)>1 else float(vValues[SelectedId[0]])
    SelectedId                  =     [a for a, x in enumerate(vNames) if x == "Sphericity"]
    Sphericity                  =     [float(vValues[id]) for id in SelectedId] if len(SelectedId)>1 else float(vValues[SelectedId[0]])
    return Volume,Sphericity
          
# Process image        
def GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing):
    global vImaris
    global numberIndex
    global SelectedChanelIndex, DAPIChannel
    global UserSetDiameters
    global ListOfContainers
    global VolumeList
    global SphericityList
    global NucleolusVolumeList
    global SphericityNucleolusList
    global IntensitySumListN     
    global IntensitySumListNW
    global NumberOfSurfacesListAll
    ListOfContainers    =   []   # To keep track of all the containers I will create in scene, so that I can remove them after saving image
    vImage              =       vImaris.GetDataSet()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        if FileIndex    == 1:
            UserSetDiameters=[]
            SelectedChanelIndex=[]
            SelectedChanelIndex, DAPIChannel= Ask_user(numberIndex)
            if len(SelectedChanelIndex)>0:
                AskUserDiameters(SelectedChanelIndex)
            if len(SelectedChanelIndex)==0 and len(UserSetDiameters)==0:
                return 1
        numberIndex     +=1
        ChromocenterSurface =   NucleusSurface  =   NucleolusSurface    =   None
        date            =    str(datetime.datetime.now()).split(" ")[0]
        date            +=    " 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        vImage                     =    vImaris.GetDataSet()
        if DoSegmentation:
            logtime('Segmentation START - image_'+str(FileIndex))
            vScene                     =    vImaris.GetSurpassScene()
            vFactory                   =    vImaris.GetFactory().CreateFactory()
            GroupOfObjects               =    vFactory.CreateDataContainer()
            GroupOfObjects.SetName('Segmented objects')
            ListOfContainers.append(GroupOfObjects)
            if DistanceOptions[1]:
                SurfNP=Segment_Surface(vImage, DAPIChannel, ParametersList[0], "Detailed Nucleus segmentation",0)
                SurfNP2=Segment_Surface(vImage, DAPIChannel, ParametersList[1], "Rough Nucleus segmentation",0)
                vImage_data             =   vImage.GetDataVolumeFloats(DAPIChannel,0)
                Mask1=Get_Mask_data(SurfNP)
                Mask2=Get_Mask_data(SurfNP2)
                resNuc=getLowIntensity(Mask1, vImage_data, Mask2)
                AddChannel(resNuc, vImage, numberIndex)
                vImage.SetChannelName(numberIndex-1,"Low DAPI intensity")
                NucR=Segment_Surface(vImage, numberIndex-1, ParametersList[2], "Rough Nucleolus segmentation",ParametersList[3])
                numberIndex+=1
                logtime('Nucleolus surface segmented - Image '+str(FileIndex))
                NucleolusSurface=SelectSurface(vScene, NucR, "Nucleolus", GroupOfObjects) #Nucleolus is segmented
            SurfP         =     Segment_Surface(vImage, DAPIChannel, ParametersList[4], "Rougher Nucleus segmentation",0)
            NucleusSurface=     SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
            logtime('Nucleus surface segmented - Image '+str(FileIndex))
            #Segment chromocenters and FISH signals
            AllChannels     =   [DAPIChannel]+SelectedChanelIndex
            ObjectList      =   ["Chromocenters"]+["FISH-Ch"+str(idx)for idx in SelectedChanelIndex]
            counter         =   0
            NumberOfSurfacesList = [] 
            SegDetail=[ParametersList[5], ParametersList[6]]+UserSetDiameters
            n=0
            for chId in AllChannels:
                ChannelName ="Channel"+str(chId)+" High intensity"
                ObjectName  =   ObjectList[counter]
                NumberOfSurfaces, TempSurf     =   SegHighIntensity(chId, SurfP, ChannelName, ObjectName, numberIndex, GroupOfObjects,SegDetail[n], SegDetail[n+1], AllChannels, Result_pathway, FileIndex)                
                if chId==DAPIChannel:
                    CCSurf      = TempSurf 
                numberIndex +=  1
                counter     +=  1
                n           +=  2
                NumberOfSurfacesList.append(NumberOfSurfaces)
            NumberOfSurfacesListAll.append(NumberOfSurfacesList)
            mask          =     GeteMasks(NucleolusSurface, CCSurf, "inverse", SurfP)
            AddChannel(mask, vImage, numberIndex)
            surfWOCCNucl            =   Segment_Surface(vImage, numberIndex-1,0.8, "SurfWithoutCC",0) #segment part of nucleus without chromocenters and nucleolus            
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface,ImmunoNames,ImmunoSignalList  = GetSegmentedObjects(DistanceOptions, FileIndex)
        #Get volume and sphericity of nucleus and nucleolus
        SurfaceVolume, Sphericity    =   GetNucParameters(NucleusSurface)
        VolumeList.append(SurfaceVolume)
        SphericityList.append(Sphericity)
        SurfaceVolume, Sphericity    =   GetNucParameters(NucleolusSurface)
        NucleolusVolumeList.append(SurfaceVolume)
        SphericityNucleolusList.append(Sphericity)
        IntensitySum     =   GetNucleusIntensity(surfWOCCNucl, AllChannels)
        IntensitySumListN.append(IntensitySum)       
        IntensitySum     =   GetNucleusIntensity(CCSurf, AllChannels)
        IntensitySumListNW.append(IntensitySum)
        if DoSegmentation:
            vPath = os.path.join(Result_pathway, vFileName+".ims")
            vImaris.FileSave(vPath, "")
    vScene.RemoveChild(GroupOfObjects)
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
  label_text    =    Messge
  option        =    OPTIONS
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
    Messge                          =       "Please select the RNA PolII immunostaining channels: \n You can choose multiple channels."
    PopUpMessage(OPTIONS, Messge)
    SelectedChanelIndex                =        [i for i, x in enumerate(User_selection) if x == 1]
    Messge="Please select the DAPI channel: \n Please only choose one channel."
    PopUpMessage(OPTIONS, Messge)
    DAPIChannel    =    [i for i, x in enumerate(User_selection) if x == 1][0]
    return SelectedChanelIndex, DAPIChannel

#==============================================================================
# Functions required to ask user set immuno signal spot diameter for all channels
#==============================================================================
def fetch(entries, root):
    global UserSetDiameters
    TempList    =  []
    for entry in entries:
        text  = entry.get()
        try:
            x = float(text)
            TempList  =   TempList + [x]
        except ValueError:
            print "You must enter a number in all the fields"
            return
    root.destroy()
    UserSetDiameters    =   TempList

def makeform(root, fields , selectedChannels):
    entries = []
    for ChannelID in selectedChannels:
        for txt in fields:
            field = txt +" for Channel "+str(ChannelID)
            row = Frame(root)
            lab = Label(row, width=30, text=field, anchor='w')
            ent = Entry(row)
            row.pack(side=TOP, fill=X, padx=5, pady=5)
            lab.pack(side=LEFT)
            ent.pack(side=RIGHT, expand=YES, fill=X)
            entries.append(ent)
    return entries

def AskUserDiameters(selectedChannels):
    fields  =    ['Diameter XY', 'Diameter Z']
    root = Tk()
    label_text="Please set the diameters for \nimmunostaining signal segmentation.\nPlease enter a number in all the fields"
    row = Frame(root)
    lab = Label(row, width=30, text=label_text, anchor='w')
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ents = makeform(root, fields, selectedChannels)
    root.bind('<Return>', (lambda event, e=ents: fetch(e, root)))
    b1 = Button(root,fg="darkgreen", text='Submit',
          command=(lambda e=ents: fetch(e, root)))
    b1.pack(side=LEFT, padx=5, pady=5)
    b2 = Button(root,fg="darkred", text='Quit', command=root.quit)
    b2.pack(side=LEFT, padx=5, pady=5)
    root.mainloop()

#==============================================================================
# Function to get parameters for this plugin from file:XTSegmentNucleus_Parmaeters.csv
#==============================================================================
def GetPluginParameters():
    currentDirectory                         =       os.getcwd()
    AllFilesInDirectory                      =       os.listdir(currentDirectory)
    ParametersList                           =       None
    if "XTSegmentNucleus_Parameters.csv" in AllFilesInDirectory:
        ParameterData                            =       pd.read_csv("XTSegmentNucleus_Parameters.csv", sep=";", header='infer',decimal='.')
        if "Value" in ParameterData.columns:
            ParametersList                           =       list(ParameterData["Value"])
        else:
            tkMessageBox.showinfo(title="Error", message="Please make sure the 'XTSegmentNucleus_Parameters.csv' file contains a column 'Value' containing the values necessary for this plugin.")
            quit()
    else:
        tkMessageBox.showinfo(title="Error", message="Please make sure there is a 'XTSegmentNucleus_Parameters.csv' in the folder containing the 'XTSegmentNucleus.py'.")
        quit()
    return ParametersList

#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTSegmentNucleus_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTSegmentNucleusTest(aImarisId):
    logging.basicConfig(level=logging.DEBUG, filename= "log[XTSegmentNucleus].log")
    try:
        #Declare global variables
        global gLasttime
        global vImaris
        global VolumeList
        global SphericityList
        global NucleolusVolumeList
        global SphericityNucleolusList
        global IntensitySumListN     
        global IntensitySumListNW
        global NumberOfSurfacesListAll
        NumberOfSurfacesListAll =   []
        VolumeList              =   []
        SphericityList              =   []
        NucleolusVolumeList     =   []
        SphericityNucleolusList =   []
        IntensitySumListN       =   []  
        IntensitySumListNW      =   []
        gLasttime               =    None
        logtime('Extension XTSegmentNucleus START')
        FileNameList            =   []
        #            Step1: Connect to Imaris
        #==============================================================================
        vImarisLib            =    ImarisLib.ImarisLib()
        # Get an imaris object with id aImarisId
        vImaris             =   vImarisLib.GetApplication(aImarisId)
        ParametersList      =   GetPluginParameters()
        # """BEGIN LOOP OVER ALL IMAGES"""
        # Open File and get filename
        if vImaris is not None :
            logtime('Connected to Imaris')
            ListOfOptions                   =       [["Batch of images", "Just one image"], ["Segmentation & Quantify RNA PolII", "Quantify RNA PolII"], ["Nucleus", "Nucleolus", "Chromocenters"]]
            ListOfMessages                  =       ["Do you wish to run the script on a batch of images or just on one image already opened?", "Do you wish to do automated segmentation?", "Do you wish to analyse RNA PolII distribution's as a function of the:"]
            UserParameterList               =       []
            for i in range(len(ListOfOptions)):
                OPTIONS                         =       ListOfOptions[i]
                Messge                          =       ListOfMessages[i]
                PopUpMessage(OPTIONS, Messge)
                UserParameterList               =       UserParameterList + [User_selection]
            BatchProcessing                    =        UserParameterList[0][0]
            DoSegmentation                    =        UserParameterList[1][0]
            DistanceOptions                 =       UserParameterList[2]
            TypeStudy           =   []   # List variable to track the type of study selected by user :  distribution as a function of the nucleus and / or chromocenters and / or nucleolus 
            TypeStudy.append("Nucleus")
            surfaceAnalised     =   [x for x,y in zip(["Chromocenters", "Nucleolus"], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
            SurfaceName         =   '-'.join(surfaceAnalised)
            TypeStudy.extend(SurfaceName)
            FileIndex                       =       1 #This  variable is used to count the number of files analysed
            if BatchProcessing  :
                #         Step2: Here the user is asked to set the path to the folder containing image to be analysed
                #==============================================================================
                root1                    =            Tk()
                Image_folder            =            tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
                root1.destroy()
                Result_pathway          =           os.path.join(Image_folder, "XTSegmentNucleus_Result")
                CreateDirectoryToSaveFiles(Result_pathway)
                AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
                logtime('Get all files')
                AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
                TotalNumberFile         =           len(AllFilesToBeProcessed)
                logtime('Start bach processing')
                if TotalNumberFile  >   0:
                    for vFileName in AllFilesToBeProcessed:
                        vFullFileName   =           os.path.join(Image_folder, vFileName)
                        vImaris.FileOpen(vFullFileName, "")
                        ImageIsEmpty    =           GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing)
                        if  not ImageIsEmpty :
                            FileIndex                   +=          1
                            FileNameList.append(vFileName)
                        else:
                            tkMessageBox.showinfo(title="Alert", message="User needs to select channels and set diameter values for spot segmentation.")
                            quit()
                else:
                    tkMessageBox.showinfo(title="Alert", message="There is no .ims or .ics file detected in the selected folder.")
                    quit()
            else    :
                TotalNumberFile         =   1
                vFileName               =   vImaris.GetCurrentFileName()
                ImageIsEmpty            =   True
                if vFileName !=""       :
                    vFileName           =   vImaris.GetCurrentFileName()
                    vImage              =   vImaris.GetDataSet()
                    if vImage is None :
                        tkMessageBox.showinfo(title="Alert", message="No image is detected. \n Please open an image and select on 'SegmentNucleus_Result' again.")
                        quit()
                    vFilePath           =   os.path.dirname(vFileName)
                    Result_pathway      =   os.path.join(vFilePath, "XTSegmentNucleus_Result")
                    CreateDirectoryToSaveFiles(Result_pathway)
                    vFileName           =   os.path.split(vFileName)[1]
                    ImageIsEmpty        =   GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing)
                if not ImageIsEmpty :
                    FileIndex                   +=          1
                    FileNameList.append(vFileName)                    
                if ImageIsEmpty:
                    tkMessageBox.showinfo(title="Alert", message="User needs to select channels and set diameter values for spot segmentation.")
                    quit()
            if FileIndex>1:
                logtime('Files saved and organised')
                PullSimilarDataIntoOneFile(FileNameList, Result_pathway)
            CreateTables(FileNameList, VolumeList, "NucleusVolume", Result_pathway, "", 0, "") 
            CreateTables(FileNameList, SphericityList, "NucleusSphericity", Result_pathway, "", 0, "")   
            CreateTables(FileNameList, SphericityNucleolusList, "SphericityNucleolusList", Result_pathway, "", 0, "")            
            CreateTables(FileNameList, IntensitySumListN, "IntensitySumListN", Result_pathway, "", 0, "")            
            CreateTables(FileNameList, IntensitySumListNW, "IntensitySumListCC", Result_pathway, "", 0, "")    
            CreateTables(FileNameList, NumberOfSurfacesListAll, "NumberOfSurfacesListAll", Result_pathway, "", 0, "")    

            logtime('XTSegmentNucleus extension done')
            print "All tasks have been completed successfully. \n Resulting files are saved in the folder XTSegmentNucleus_Result"
            raw_input("Press Enter to terminate.")
        else:
            tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
            logtime('Extension XTSegmentNucleus END')
    except:
        logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================