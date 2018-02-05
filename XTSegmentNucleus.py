# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
# Note: This script is calibrated for 3D images of plant nuclei obtained using Leica TCS SP8
# Creator: Mariamawit S. Ashenafi
# Published on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="SegmentNucleus" icon="Python" tooltip="SegmentNucleus">
#         <Command>PythonXT::XTSegmentNucleus(%i)</Command>
#       </Item>
#      </Menu>
#    </CustomTools>

import numpy as np
import pandas as pd
import logging
import tkMessageBox
import os
from Tkinter import *
import tkFileDialog
import time
import datetime
import ImarisLib

#==============================================================================
# Start of extension
#==============================================================================
                
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
#Function to export output in .csv files
   
#Funtion to get mean intensity inside nucleus surface for all channels selected by the user
def GetNucleusIntensity(surf, Dyes):
    AllChannel                =   [DAPIChannel]   +   Dyes
    vAllStatistics      =     surf.GetStatistics()
    vNames              =     vAllStatistics.mNames
    vValues             =     vAllStatistics.mValues
    IntenityList        =     [float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Sum"]
    IntensitySum       =   [IntenityList[y] for y in AllChannel]
    return IntensitySum

# Function to get channel intensities in all channels in nucleus shells:
def get_intensity(Dyes, Surf, NumberSurf):
    vAllStatistics               =     Surf.GetStatistics()
    vNames                       =     vAllStatistics.mNames
    vValues                      =     vAllStatistics.mValues
    #Intensity properties for each spots
    label                        =    "Intensity Sum"
    Index                        =    [i for i, x in enumerate(vNames) if x == label]
    result                       =    [vValues[e] for e in Index]
    IntensityInAllChannel        =    [result[(vChannelI*NumberSurf):(vChannelI*NumberSurf+NumberSurf)] for vChannelI in Dyes]
    return IntensityInAllChannel
      
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
    vImage            =    vImaris.GetDataSet()
    vImageSizeX     =    vImage.GetSizeX()
    vImageSizeY     =    vImage.GetSizeY()
    vImageSizeZ     =    vImage.GetSizeZ()
    vExtentMinX        =    vImage.GetExtendMinX()
    vExtentMinY        =    vImage.GetExtendMinY()
    vExtentMinZ        =    vImage.GetExtendMinZ()
    vExtentMaxX        =    vImage.GetExtendMaxX()
    vExtentMaxY        =    vImage.GetExtendMaxY()
    vExtentMaxZ        =    vImage.GetExtendMaxZ()
    mask_min        =    [vExtentMinX, vExtentMinY, vExtentMinZ]
    mask_max        =    [vExtentMaxX, vExtentMaxY, vExtentMaxZ]
    mask_size        =    [vImageSizeX, vImageSizeY, vImageSizeZ]
    mask_time        =    0
    mask            =    Surf.GetMask(mask_min[0], mask_min[1], mask_min[2], mask_max[0], mask_max[1], mask_max[2],mask_size[0],mask_size[1], mask_size[2], mask_time)
    mask_values        =    mask.GetDataVolumeFloats(0,0)
    return mask_values

#This function, inside a given surface, selectes voxels containing 50% of the highest intensity, sets there value to the maximum intensity
def SelectVoxels(mask_values, vImage_data, Type):
    maxValue=getMax(vImage_data)
    Result                  =   []
    for z in range(len(vImage_data)):
        maskY               =   []
        for y in range(len(vImage_data[0])):
            mask_surface    =    mask_values[z][y]
            mask            =    vImage_data[z][y]
            if Type=="HighIntensity":
                mask            =   [round(item/maxValue, 1) for item in mask]
                mask            =    [maxValue if item>0.4 else item for item in mask]
            mask            =    [x*y for x,y in zip(mask, mask_surface)]
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
    vSurpassScene         =     vImaris.GetSurpassScene()
    vFactory              =     vImaris.GetFactory()
    number_scene_instance =     range(vSurpassScene.GetNumberOfChildren())
    volume_notFound       =     True
    i                     =     0
    while i <= number_scene_instance and volume_notFound:
        selection     =    vSurpassScene.GetChild(i)
        vVolume     =    vFactory.ToVolume(selection)
        if not vVolume is None:
            volume_notFound=False
        i=i+1
    vAllStatistics         =     vVolume.GetStatistics()
    vNames                   =     vAllStatistics.mNames
    vValues                =     vAllStatistics.mValues
    Max_intensity_index    =     [a for a, x in enumerate(vNames) if x == "Data Intensity Max"]
    Max_intensity        =    [vValues[x] for x in Max_intensity_index]
    return Max_intensity


#Segments a channel into surface and returns the surface
def Segment_Surface(vImage, ch, vSFW, name,vLCF, vSFS):
    vROI = None
    vATA = 1
    vATM = 0
#    vSFS = ''
    vDCI = ch
    vSurface2 = vImaris.GetImageProcessing().DetectSurfaces(vImage, vROI, vDCI,vSFW, vLCF, vATA, vATM, vSFS)
    vSurface2.SetName(name)
    return vSurface2

#This function, inside a given surface, selectes voxels containing 50% of the lowest intensity, sets there value to the maximum intensity, and sets the value of the other voxels to 0
def getLowIntensity(mask_values, vImage_data, mask_values2):
    maxValue                =   getMax(vImage_data)
    Result                  =   []
    for z in range(len(vImage_data)):
        mask_y              =   []
        for y in range(len(vImage_data[0])):
            mask_surface    =    mask_values[z][y]
            mask_surface2    =    mask_values2[z][y]
            mask            =    vImage_data[z][y]
            # mask=[1 for item==0 else item for item in mask]
            mask            =    [round(item/maxValue, 1) for item in mask]
            mask            =   [0 if item>=0.5 else item for item in mask]
            mask_surface    =   [maxValue if item==0 else 0 for item in mask_surface]
            mask            =    [i*j*k for i,j,k in zip(mask, mask_surface, mask_surface2)]
            mask_y.append(mask)
        Result.append(mask_y)
    return Result

#For a surface containing several surface, this function selectes the most dense surface and creates a surface containing only this selected surface
def SelectSurface(vscene, surf, surfType, groupContainer):
    vTimeIndex          =   0
    vAllStatistics         =     surf.GetStatistics()
    vNames                   =     vAllStatistics.mNames
    vValues                =     vAllStatistics.mValues
    volID                =     [a for a, x in enumerate(vNames) if x == "Volume"]
    vol                    =    [vValues[x] for x in volID]
    if surfType         ==  "Nucleolus":
        areaID            =     [a for a, x in enumerate(vNames) if x == "Area"]
        area            =    [vValues[x] for x in areaID]
        volumeToAreaRatio         =   [x/y for x,y in zip(vol, area)]
        SelectedID      =   max(xrange(len(volumeToAreaRatio)), key=volumeToAreaRatio.__getitem__)
    else:
        SelectedID      =   max(xrange(len(vol)), key=vol.__getitem__)
    verts=surf.GetVertices(SelectedID)
    vNormals            =    surf.GetNormals(SelectedID)
    faces=surf.GetTriangles(SelectedID)
    vNucleiSurface        =    vImaris.GetFactory().CreateSurfaces()
    vNucleiSurface.AddSurface(verts, faces, vNormals, vTimeIndex)
    vNucleiSurface.SetName(surfType)
    groupContainer.AddChild(vNucleiSurface, -1)
    vscene.AddChild(groupContainer, -1)
    return vNucleiSurface

#This function first creates a channel by selecting voxel contianing 50% of highest intensity values
def SegHighIntensity(chIndex, surf, ChannelName, ObjectName, numberIndex, GroupContainer, SmoothDetail, ObjectDiam, Dyes, Result_pathway, FileIndex):
    vImage                  =   vImaris.GetDataSet()
    vScene                  =   vImaris.GetSurpassScene()
    vImage_data             =   vImage.GetDataVolumeFloats(chIndex,0)
    MaskSurface             =   Get_Mask_data(surf)
#    result                  =   SelectVoxels(MaskSurface, vImage_data, "HighIntensity")
    result                  =   SelectVoxels(MaskSurface, vImage_data, "")
    AddChannel(result, vImage, numberIndex) #CC channel
    vImage.SetChannelName(numberIndex-1,ChannelName)
    Surface                 =   Segment_Surface(vImage, numberIndex-1, SmoothDetail, ObjectName,ObjectDiam, "")
    GroupContainer.AddChild(Surface, -1)
    vScene.AddChild(GroupContainer, -1)
    NumberOfSurfaces        =    Surface.GetNumberOfSurfaces()
    IntenistyList           =   get_intensity(Dyes, Surface, NumberOfSurfaces)
    positionOfAllItems      =   [Surface.GetCenterOfMass(surfId)[0] for surfId in range(NumberOfSurfaces)]
    SurfVolume              =       GetNucParameters(Surface, "Volume")
#    if type(SurfVolume) is float:
#        SurfVolume          =   [SurfVolume]
    CreateTables(None, IntenistyList, "SurfaceIntensity_Ch", Result_pathway, chIndex, 1, FileIndex)  # export intensities in selected channels for each spot
    CreateTables(None, positionOfAllItems, "SurfacePositions_Ch", Result_pathway, chIndex, 0, FileIndex)
    CreateTables(None, SurfVolume, "SurfaceVolume_Ch", Result_pathway, chIndex, 0, FileIndex)
    return NumberOfSurfaces, Surface

def GeteMasks(NucleolusSurfaace, CCSurf, Type, NuucleusSurface):
    NucMask             =       Get_Mask_data(NuucleusSurface) if NuucleusSurface is not None else None
    mask_values         =       Get_Mask_data(NucleolusSurfaace) if NucleolusSurfaace is not None else None
    mask_values1         =       Get_Mask_data(CCSurf) if CCSurf is not None else None
    if Type=="inverse":
        MaskX           =   []
        for x in range(len(mask_values)):
            MaskY       =   []
            for y in range(len(mask_values[x])):
                MaskZ   =   mask_values[x][y]
                MaskCC  =   mask_values1[x][y]
                MaskNuc =   NucMask[x][y]
                MaskZ   =   [b if a==0 and c==0 else 0 for a,b, c in zip(MaskZ, MaskNuc, MaskCC)]
                MaskY.append(MaskZ)
            MaskX.append(MaskY)
        mask_values   =   MaskX
    return mask_values


#This function returns the volume value of surface
def GetNucParameters(surf, LabelType):
    vAllStatistics              =     surf.GetStatistics()
    vNames                      =     vAllStatistics.mNames
    vValues                     =     vAllStatistics.mValues
    SelectedId                  =     [a for a, x in enumerate(vNames) if x == LabelType] 
    Result                      =     [float(vValues[id]) for id in SelectedId] 
    return Result
          
# Process image        
def GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName):
    global vImaris
    global numberIndex
    global SelectedChanelIndex, DAPIChannel
    global UserSetDiameters
    global ListOfContainers
    global IntensitySumListN     
    global IntensitySumListNW
    global NumberOfSurfacesListAll
    global OutputList
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
            SurfP         =     Segment_Surface(vImage, DAPIChannel, ParametersList[4], "Rougher Nucleus segmentation",0, "")
            NucleusSurface=     SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
            if DistanceOptions[1]:
                SurfNP=Segment_Surface(vImage, DAPIChannel, ParametersList[0], "Detailed Nucleus segmentation",0, "")
                SurfNP2=Segment_Surface(vImage, DAPIChannel, ParametersList[1], "Rough Nucleus segmentation",0, "")
                vImage_data             =   vImage.GetDataVolumeFloats(DAPIChannel,0)
                Mask1=Get_Mask_data(SurfNP)
                Mask2=Get_Mask_data(SurfNP2)
                resNuc=getLowIntensity(Mask1, vImage_data, Mask2)
                AddChannel(resNuc, vImage, numberIndex)
                vImage.SetChannelName(numberIndex-1,"Low DAPI intensity")
                NucleusArea     =    GetNucParameters(SurfP, "Area")[0]/2.0
#                SurfThres   =   '"Area" between 1.00 um^2 and '+str(NucleusArea)+' um^2'
                SurfThres   =""
                NucR=Segment_Surface(vImage, numberIndex-1, ParametersList[2], "Rough Nucleolus segmentation",ParametersList[3], SurfThres)
                NucleolusSurface=SelectSurface(vScene, NucR, "Nucleolus", GroupOfObjects) #Nucleolus is segmented
#                NucleolusSurface.SetName("Nucleolus")
#                GroupOfObjects.AddChild(NucleolusSurface, -1)
                numberIndex+=1
                logtime('Nucleolus surface segmented - Image '+str(FileIndex))
            logtime('Nucleus surface segmented - Image '+str(FileIndex))
            #Segment chromocenters and FISH signals
            AllChannels     =   [DAPIChannel]+SelectedChanelIndex
            ObjectList      =   ["Chromocenters"]+["FISH-Ch"+str(idx)for idx in SelectedChanelIndex]
            NumberOfSurfacesList = [] 
            SegDetail=[ParametersList[5], ParametersList[6]]+UserSetDiameters
            n=0
            CreatedFISHChannel =    []
            for chId, ObjectName in zip(AllChannels, ObjectList):
                ChannelName ="Channel"+str(chId)+" High intensity"
                NumberOfSurfaces, TempSurf     =   SegHighIntensity(chId, SurfP, ChannelName, ObjectName, numberIndex, GroupOfObjects,SegDetail[n], SegDetail[n+1], AllChannels, Result_pathway, FileIndex)                
                if chId==DAPIChannel:
                    CCSurf      = TempSurf 
                else:
                    CreatedFISHChannel.append(numberIndex)
                numberIndex +=  1
                n           +=  2
                NumberOfSurfacesList.append(NumberOfSurfaces)
            NumberOfSurfacesListAll.append(NumberOfSurfacesList)
            mask          =     GeteMasks(NucleolusSurface, CCSurf, "inverse", SurfP)
            AddChannel(mask, vImage, numberIndex)
            surfWOCCNucl            =   Segment_Surface(vImage, numberIndex-1,0.8, "SurfWithoutCC",0, "") #segment part of nucleus without chromocenters and nucleolus            
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface,ImmunoNames,ImmunoSignalList  = GetSegmentedObjects(DistanceOptions, FileIndex)
        #Get volume and sphericity of nucleus and nucleolus
        SurfList=[NucleusSurface, NucleusSurface, NucleolusSurface, NucleolusSurface]
        TypeOfLabel=["Volume", "Sphericity"]*2         
        for lab, surf, Outp in zip(TypeOfLabel, SurfList, OutputList):
            OutPutResult    =   GetNucParameters(surf, lab)
            Outp.append(OutPutResult)
        IntensitySum     =   GetNucleusIntensity(surfWOCCNucl, CreatedFISHChannel)
        IntensitySumListN.append(IntensitySum)       
        IntensitySum     =   GetNucleusIntensity(CCSurf, CreatedFISHChannel)
        IntensitySumListNW.append(IntensitySum)
        if DoSegmentation:
            vPath = os.path.join(Result_pathway, vFileName+".ims")
            vImaris.FileSave(vPath, "")
            os.remove(vFullFileName)
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
def XTSegmentNucleus(aImarisId):
    logging.basicConfig(level=logging.DEBUG, filename= "log[XTSegmentNucleus].log")
    try:
        #Declare global variables
        global gLasttime
        global vImaris
        global IntensitySumListN     
        global IntensitySumListNW
        global NumberOfSurfacesListAll
        global OutputList
        OutputList  =   [[]]*4
        NumberOfSurfacesListAll =   []
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
                        ImageIsEmpty    =           GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName)
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
            OutputListb=[IntensitySumListN,  IntensitySumListNW, NumberOfSurfacesListAll]
            for i in OutputListb:
                OutputList.append(i)
            OutPutLabel=["NucleusVolume", "NucleusSphericity", "NucleolusVolume", "SphericityNucleolusList", "IntensitySumListN", "IntensitySumListCC", "NumberOfSurfacesListAll"]
            for tab,  lab in zip(OutputList, OutPutLabel):
                CreateTables(FileNameList, tab, lab, Result_pathway, "", 0, "") 
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