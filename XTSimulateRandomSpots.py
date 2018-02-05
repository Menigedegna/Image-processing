# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
# 	 Create masks using these surfaces (Mask1: nucleus, Mask2: nucleus without chromocenters and nucleolus),       
#   Simulate random spots in side these masks by setting the total number of spots,
#   Get DAPI intensity in each spots,
#   Exports result into .csv tables and save .ims file containing the surfaces, spots and masks created
# Note: This script is calibrated for 3D images of plant nuclei conternstaining with DPAI, obtained using Leica TCS SP8
# Creator: Mariamawit S. Ashenafi, Célia Baroux, UZH
# Published on 08.06.2017
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTSimulateRandomSpots" icon="Python" tooltip="XTSimulateRandomSpots">
#         <Command>PythonXT::XTSimulateRandomSpots(%i)</Command>
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
import XTSegmentNuclei as SN


#==============================================================================
# Start of extension
#==============================================================================
# Function to generate random spots in nucleus
def simulate3DPointsInSurface(vFactory, groupContainer, aRadius, chID, number_spots, surfMask, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMaxX, vExtentMaxY, vExtentMaxZ):
    """Choose random voxel IDs from mask created"""
    NumberOfSimulation              =   min(len(surfMask), number_spots)
    RandomPointVoxelId        =   np.random.choice(range(len(surfMask)), NumberOfSimulation)
    RandomPointVoxel          =   [surfMask[x] for x in RandomPointVoxelId]
    """Concert voxel ID into X,Y,Z position coordiantes"""
    RandomPosition      =   ConvertCoordianteIntoVoxelId("Rev", RandomPointVoxel, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMaxX, vExtentMaxY, vExtentMaxZ)
#    SpotName            =   str(NumberOfSimulation)+"random spots created"
#    vSpot               =   CreateSpots(vFactory, RandomPosition, SpotName, groupContainer, aRadius)
    return RandomPosition
 
def GetRandomIntensity(FlattenedSurfMask, FlattenedChannelData, NumbRandomSpots):
    RandomVoxelIntensityAllChannels=[]
    vImage                  =   vImaris.GetDataSet()
    VoxelsInsideSurface     =   [x for x,y in zip(range(len(FlattenedSurfMask)), FlattenedSurfMask) if y==1] #Collect ids of voxels inside mask 
    NumberOfSimulation      =   min(len(VoxelsInsideSurface), NumbRandomSpots)
    RandomVoxelList         =   np.random.choice(VoxelsInsideSurface, NumberOfSimulation)
    for DataChannelX in FlattenedChannelData:
        RandomDataChannelX  =   [DataChannelX[x] for x in RandomVoxelList]
        RandomVoxelIntensityAllChannels.append(RandomDataChannelX)
    return RandomVoxelIntensityAllChannels
    
# Function to create spots
def CreateSpots(vFactory, aPositionsXYZ, SpotName, groupContainer, aRadius):
    vSpot			        =	vFactory.CreateSpots()
    aIndicesT				=	[0.0]*len(aPositionsXYZ)
    aRadii					=	[aRadius]*len(aPositionsXYZ)
    vSpot.Set(aPositionsXYZ,aIndicesT,aRadii)
    vSpot.SetName(SpotName)
    groupContainer.AddChild(vSpot, -1)
    return vSpot

#Funtion to get mean intensity inside nucleus surface for all channels selected by the user
def GetNucleusIntensity(surf):
    Dyes                =   [DAPIChannel]   +   SelectedChanelIndex
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    IntenityList        = 	[float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Mean"]
    IntenityListSel     =   [IntenityList[y] for y in Dyes]
    return IntenityListSel

# Function to get channel intensities in all channels in nucleus shells:
def get_intensity(Dyes, vSpots, Number_spot, MeanIntensitySurfaceList):
    vAllStatistics 				= 	vSpots.GetStatistics()
    vNames	      				= 	vAllStatistics.mNames
    vValues      			    = 	vAllStatistics.mValues
    #Intensity properties for each spots
    label					    =	"Intensity Sum"
    Index			          	=	[i for i, x in enumerate(vNames) if x == label]
    result			            =	[vValues[e] for e in Index]
    IntensityInAllChannel       =   [result[(vChannelI*Number_spot):(vChannelI*Number_spot+Number_spot)] for vChannelI in Dyes]
    NormalisedIntensityList     =   []
    for ch in range(len(MeanIntensitySurfaceList)):
        NormalisedIntensity         =   [round(float(z)/float(MeanIntensitySurfaceList[ch]), 3) if float(MeanIntensitySurfaceList[ch])>0 else 0 for z in IntensityInAllChannel[ch]]
        NormalisedIntensityList.append(NormalisedIntensity)
    return NormalisedIntensityList

# Function to convert spots coordinates into voxel index
def ConvertCoordianteIntoVoxelId(Type, SpotPosition, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMaxX, vExtentMaxY, vExtentMaxZ):
    Result      =   []
    if Type=="Fwd":
        for spot in SpotPosition:
            VoxelID     =   [int(round(x*y, 0)) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
            Result.append(VoxelID)
    else:
        for spot in SpotPosition:
            Coordinates     =   [round(x*y, 3) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
            Result.append(Coordinates)
    return Result


#This function returns the volume value of surface
def GetVolume(surf):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    VolumeId	           = 	[a for a, x in enumerate(vNames) if x == "Volume"]
    Volume		           =	float(vValues[VolumeId[0]])
    return Volume

def GeteMasks(ListSurface, Type, NuucleusSurface):
    ListSurfaceMasks   =   []
    NucMask         =   SN.Get_Mask_data(NuucleusSurface) if NuucleusSurface is not None else None
    for SurfaceInstance in ListSurface:
        mask_values                 =   SN.Get_Mask_data(SurfaceInstance)
        if Type=="inverse":
            MaskX       =   []
            for x in range(len(mask_values)):
                MaskY       =   []
                for y in range(len(mask_values[x])):
                    MaskZ   =   mask_values[x][y]
                    MaskNuc =   NucMask[x][y]
                    MaskZ   =   [b if a==0 else 0 for a,b in zip(MaskZ, MaskNuc)]
                    MaskY.append(MaskZ)
                MaskX.append(MaskY)
            mask_values=MaskX
        ListSurfaceMasks.append(mask_values)
    return ListSurfaceMasks

def CreateShells(vFactory, SurfaceList, number_of_shell, SurfaceName, MaxFactor):
    global ListOfContainers
    vSurpassScene 		=   vImaris.GetSurpassScene()
    ListOfShells        =   []
    GroupOfShell        =	vFactory.CreateDataContainer()
    GroupOfShell.SetName(SurfaceName + ' shells')
    ListOfContainers.append(GroupOfShell)
    vTimeIndex			=	0
    number_of_shell     =   int(number_of_shell)
#    FactorIncrease      =   0
#    if SurfaceName  !=  "Nucleus":
#        Stepper           =   (float(MaxFactor)-1.0)/float(number_of_shell) # increasing factor
    for shell in range(number_of_shell, 0, -1): # for each shell to be created
        vNewVerticesTot             =    vTrianglesTot   =   vNormalsTot =   []
        NumberOfVerticesPerSurf     =    NumberOfTrianglessPerSurf  =   vTimeIndexTot   =   []
        NewShellSurface		        =	 vImaris.GetFactory().CreateSurfaces()
        for surf in SurfaceList: 
            NumberOfSurfaces            =    surf.GetNumberOfSurfaces()
            for surfID in range(NumberOfSurfaces): # for all surfaces included in chromocenter surface
                vNormals			=	surf.GetNormals(surfID)
                vTriangles			=	surf.GetTriangles(surfID)
                vVertices			=	surf.GetVertices(surfID)
                vCenterOfMass		=	surf.GetCenterOfMass(surfID)[0]
                if SurfaceName  ==  "Nucleus":
                    vNewVertices        =   [[(x[0]-vCenterOfMass[0])*(float(shell)/number_of_shell)+vCenterOfMass[0], (x[1]-vCenterOfMass[1])*(float(shell)/number_of_shell)+vCenterOfMass[1], (x[2]-vCenterOfMass[2])*(float(shell)/number_of_shell)+vCenterOfMass[2]] for x in vVertices]
                else:
#                    FactorIncrease      +=  Stepper
                    vNewVertices        =   [[(x[0]-vCenterOfMass[0])*((number_of_shell+float(shell))/number_of_shell)+vCenterOfMass[0],(x[1]-vCenterOfMass[1])*((number_of_shell+float(shell))/number_of_shell)+vCenterOfMass[1], (x[2]-vCenterOfMass[2])*((number_of_shell+float(shell))/number_of_shell)+vCenterOfMass[2]] for x in vVertices]
                vNewVerticesTot         =   vNewVerticesTot +   vNewVertices
                vTrianglesTot           =   vTrianglesTot   +   vTriangles
                vNormalsTot             =   vNormalsTot     +   vNormals
                NumberOfVerticesPerSurf     =    NumberOfVerticesPerSurf    +   [len(vNewVertices)]
                NumberOfTrianglessPerSurf   =    NumberOfTrianglessPerSurf    +   [len(vTriangles)]
            vTimeIndexTot           =    vTimeIndexTot  +   [vTimeIndex]*NumberOfSurfaces
        NewShellSurface.AddSurfacesList(vNewVerticesTot, NumberOfVerticesPerSurf,  vTrianglesTot, NumberOfTrianglessPerSurf,  vNormalsTot, vTimeIndexTot)
        NewShellSurface.SetName(SurfaceName+" shell_"+str(shell))
        GroupOfShell.AddChild(NewShellSurface, -1)
        ListOfShells        =   ListOfShells+[NewShellSurface]
    vSurpassScene.AddChild(GroupOfShell, -1)
    return ListOfShells

#==============================================================================
# These functions are required to do spot distribution analysis :
# Create shells for nucleus, for chromocenters and or for nucleolus
# Get masks for each shell
# Use masks to identify voxels inside shells
# Randomly select voxels and get X,Y,Z coordinates and intensity in selected channels
#==============================================================================
def SimulateRandomSpots(TypeStudy, numberIndex, ListParameters, ChromocenterSurface, NucleusSurface, NucleolusSurface, NucleusSurfaceVolume, FileIndex, Result_pathway, DistanceOptions):
    global ListOfContainers 
    Dyes                        =   [DAPIChannel]   +   SelectedChanelIndex # Indexes of intensities I save for each spot
    NumberSelectedChannel       =   len(SelectedChanelIndex)
    vImage			         =   vImaris.GetDataSet()
    vFactory		         =   vImaris.GetFactory()
    vSurpassScene 		   =   vImaris.GetSurpassScene()
    NumberOfShell               =   1
    MeanIntensitySurfaceList    =   GetNucleusIntensity(NucleusSurface) # get intensities of all channels inside nucleus surface
    ListMask                    =   [] # Track mask of shells created
    vExtentMinX		        =	vImage.GetExtendMinX()
    vExtentMinY		        =	vImage.GetExtendMinY()
    vExtentMinZ		        =	vImage.GetExtendMinZ()
    vExtentMaxX		        =	vImage.GetExtendMaxX()
    vExtentMaxY		        =	vImage.GetExtendMaxY()
    vExtentMaxZ		        =	vImage.GetExtendMaxZ()
    SixeX                   =   vExtentMaxX -   vExtentMinX
    SixeY                   =   vExtentMaxY -   vExtentMinY
    SixeZ                   =   vExtentMaxZ -   vExtentMinZ
    vImageSizeX 	        =	vImage.GetSizeX() 
    vImageSizeY 	        =	vImage.GetSizeY()
    vImageSizeZ 	        =	vImage.GetSizeZ()
    VoxelSizeX              =   SixeX/vImageSizeX
    VoxelSizeY              =   SixeY/vImageSizeY
    VoxelSizeZ              =   SixeZ/vImageSizeZ
    logtime('Simulate random spots START')
    Label                  =    ["Obs", "Sim"] #2 type of spots are created in scene: randome simulated spots and segmented spots of selected channels
#==============================================================================
#   		Create  mask of nucleus surface, mask0
#==============================================================================
    if DistanceOptions[0] and NucleusSurface is not None:
        ShellMaskList       =       GeteMasks([NucleusSurface], "", None)[0]
        ListMask.append(ShellMaskList)
        logtime("Shells are created for Nucleus")
#==============================================================================
#   		Create mask of nucleus surface excluding chromocenters and nucleolus, mask1 
#==============================================================================
    SurfaceList             =   [x for x,y in zip([ChromocenterSurface, NucleolusSurface], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
    if len(SurfaceList)>0:
        surfaceAnalised     =   [x for x,y in zip(["Chromocenters", "Nucleolus"], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
        SurfaceName         =   '-'.join(surfaceAnalised)
        ListOfShells        =   CreateShells(vFactory, SurfaceList, NumberOfShell, "Nucleus", ListParameters[10])
        ShellMaskList       =   GeteMasks(ListOfShells, "inverse", NucleusSurface)[0]
        ListMask.append(ShellMaskList)
        logtime("Shells are created for"+SurfaceName)
    GroupOfRandomPoints                 =   vFactory.CreateDataContainer() # Container for simulated spots
    GroupOfRandomPoints.SetName("Random spots")
    ListOfContainers.append(GroupOfRandomPoints)
    NumberOfSimulation      =   3
#==============================================================================
#     Get flattened list of intensity inside each voxel for all channel and for surface mask
#==============================================================================
    FlattenedListChannel          =  []  
    for channelId in Dyes:
        ChannelData               =   vImage.GetDataVolumeFloats(channelId,0)
        FlattenedChannelData      =  [x for sublistZ in ChannelData  for sublistY in sublistZ for x in sublistY]
        FlattenedListChannel.append(FlattenedChannelData)
#==============================================================================
#         Get random spots for each mask 
#==============================================================================
    for mk in range(len(ListMask)):
        ShellMaskList       =   ListMask[mk]
        ListSelectedVoxels    =   [[x,y,z] for x in range(vImageSizeX) for y in range(vImageSizeY) for z in range(vImageSizeZ) if ShellMaskList[x][y][z]]
        FlattenedSurfMask       =   [x for sublistZ in ShellMaskList for sublistY in sublistZ for x in sublistY]
        for NumbSpots in [1000, 10000,100000,200000,5000,10000]:
            for sim in range(NumberOfSimulation):
                #Create random spots
                RandomSpot          =   simulate3DPointsInSurface(vFactory, GroupOfRandomPoints,  0.025, 0, NumbSpots, ListSelectedVoxels,  VoxelSizeX, VoxelSizeY, VoxelSizeZ, vExtentMaxX, vExtentMaxY, vExtentMaxZ)
                SaveTable         =       pd.DataFrame(RandomSpot).T
                vPathToSaveTables = os.path.join(Result_pathway, "Position_SP"+str(0)+"_ST"+str(mk)+"_Sml"+str(sim)+"_Sim_"+str(NumbSpots))
                SaveTable.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
                logtime("Got position for "+str(NumbSpots)+" simulated spots")
                #Select random voxels and get intensities in all channels
                SpotIntensity     =   GetRandomIntensity(FlattenedSurfMask, FlattenedListChannel, NumbSpots)
                SaveTable         =   pd.DataFrame(SpotIntensity).T
                SaveTable.columns =   ["DAPI"]+["Channel_"+str(chID) for chID in SelectedChanelIndex]
                vPathToSaveTables = os.path.join(Result_pathway, "Intensity_SP"+str(0)+"_ST"+str(mk)+"_Sml"+str(sim)+"_Sim_"+str(NumbSpots))
                SaveTable.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
                logtime("Got intensities for "+str(NumbSpots)+" selected voxels")
    vSurpassScene.AddChild(GroupOfRandomPoints, -1)
    # 		Save dataframes containing results
# ==============================================================================
    logtime('Simulate random spots END')

#==============================================================================
# Function required to get segmented objects
#==============================================================================
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
        ImmunoSpotNames        =    []
        ImmunoSpotList         =    []
        while i <= numberSceneInstance :
            selection 	          =    vContainer.GetChild(i)
            vObject               =    vFactory.ToSurfaces(selection)
            if vObject is not None:
                if vObject.GetName() == "Nucleus" and DistanceOptions[0]        :   NucleusSurface      = vObject
                if vObject.GetName() == "Nucleolus" and DistanceOptions[1]      :   NucleolusSurface    = vObject
                if vObject.GetName() == "Chromocenters" and DistanceOptions[2]  :   ChromocenterSurface = vObject
            vObject               =    vFactory.ToSpots(selection)
            if vObject is not None:
                PositionTable =     vObject.GetPositionsXYZ()
                if "Immuno" in vObject.GetName() and PositionTable is not None:
                    ImmunoSpotNames     = ImmunoSpotNames   +    [vObject.GetName()]
                    ImmunoSpotList      = ImmunoSpotList    +    [vObject]
            i+=1
        logtime('Object detection END - image_'+str(FileIndex))
    return NucleusSurface, NucleolusSurface, ChromocenterSurface, ImmunoSpotNames, ImmunoSpotList


#==============================================================================
# This function: removes all objects created in scene
#==============================================================================
def RemoveObjectsCreated(vScene, ListOfContainers):
    for i in ListOfContainers:
        vScene.RemoveChild(i)

#==============================================================================
# This function:
# Segments chanenls into surface (nucleus, nucleolus and chromocenters) or spots (RNA PolII foci),
# Selects random voxels to get spot position and intensity
# And saves position and intensity in csv file
#==============================================================================
def GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName):
    global vImaris
    global numberIndex
    global SelectedChanelIndex, DAPIChannel
    global UserSetDiameters
    global ListOfContainers
    ListOfContainers    =   []   # Keep track of all the containers, they will be removed before closing this image
    ChromocenterSurface =   NucleusSurface  =   NucleolusSurface    =   None
    vImage              =   	vImaris.GetDataSet()
    vScene                     =    vImaris.GetSurpassScene()
    vFactory		           =	vImaris.GetFactory().CreateFactory()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        if FileIndex    == 1:
            SelectedChanelIndex= Ask_user(numberIndex)
            if len(SelectedChanelIndex)==0:
                tkMessageBox.showinfo(title="Error", message="Please run the plugin again and make sure you have selected channels where intensities should be measured.")
                quit()
        numberIndex     +=1
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
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
            logtime('Segmentation START')
            IsImageCorrect, NucleusSurface,ChromocenterSurface,  NucleolusSurface, DAPIChannel, GroupOfObjects=SN.SegmentAndGetFeatures(vImage, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, True)
            ListOfContainers.append(GroupOfObjects)
            logtime('Segmentation END')
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface,ImmunoNames,ImmunoSignalList  = GetSegmentedObjects(DistanceOptions, FileIndex)
#==============================================================================
#             Get random position inside created spots, and get intensities in selected channels
#==============================================================================
        if NucleusSurface.GetNumberOfSurfaces()>0:
            NucleusSurfaceVolume    =   GetVolume(NucleusSurface)
            SimulateRandomSpots(TypeStudy, numberIndex, ParametersList, ChromocenterSurface, NucleusSurface, NucleolusSurface, NucleusSurfaceVolume, FileIndex, Result_pathway, DistanceOptions)
#        if DoSegmentation:
#            vPath = os.path.join(Result_pathway, vFileName+".ims")
#            vImaris.FileSave(vPath, "")
    if len(ListOfContainers)>0 and BatchProcessing:
        RemoveObjectsCreated(vScene, ListOfContainers)
#    os.remove(vFullFileName)
    logtime('Image processing END')
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
    Messge                          =       "Please select the RNA PolII immunostaining channels: \n You can choose multiple channels."
    PopUpMessage(OPTIONS, Messge)
    SelectedChanelIndex	            =	    [i for i, x in enumerate(User_selection) if x == 1]
    return SelectedChanelIndex

#==============================================================================
# Function to get parameters for this plugin from file:XTSimulateRandomSpots_Parmaeters.csv
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
            tkMessageBox.showinfo(title="Error", message="Please make sure the 'XTSimulateRandomSpots_Parameters.csv' file contains a column 'Value' containing the values necessary for this plugin.")
            quit()
    else:
        tkMessageBox.showinfo(title="Error", message="Please make sure there is a 'XTSimulateRandomSpots_Parameters.csv' in the folder containing the 'XTSimulateRandomSpots.py'.")
        quit()
    return ParametersList

#==============================================================================
# Function required to create folder to save the plugins results
#==============================================================================
#Function to create a folder under the same directory as the images to save files that are produced
def CreateDirectoryToSaveFiles(Result_pathway):
    if os.path.exists(Result_pathway):
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTSimulateRandomSpots_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTSimulateRandomSpots(aImarisId):
	logging.basicConfig(level=logging.DEBUG, filename= "log[XTSimulateRandomSpots].log")
	try:
		#Declare global variables
		global gLasttime
		global vImaris
		global SelectedChanelIndex
		gLasttime               =    None
		logtime('Extension XTSimulateRandomSpots START')
		FileNameList            =   []
		#            Step1: Connect to Imaris
		#==============================================================================
		vImarisLib			=	ImarisLib.ImarisLib()
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
			ListOfOptions                   =       [["Batch of images", "Just one image"], ["Segmentation & Quantify RNA PolII", "Quantify RNA PolII"], ["Nucleus", "Nucleolus", "Chromocenters"]]
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
			TypeStudy           =   []   # List variable to track the type of study selected by user :  distribution as a function of the nucleus and / or chromocenters and / or nucleolus 
			TypeStudy.append("Nucleus")
			surfaceAnalised     =   [x for x,y in zip(["Chromocenters", "Nucleolus"], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
			SurfaceName         =   '-'.join(surfaceAnalised)
			TypeStudy.extend(SurfaceName)
			FileIndex                       =       1 #This  variable is used to count the number of files analysed
			if BatchProcessing  :
				# 		Step2: Here the user is asked to set the path to the folder containing image to be analysed
				#==============================================================================
				root1				    =	        Tk()
				Image_folder			=	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
				root1.destroy()
#				Result_pathway          =           os.path.join(Image_folder, "XTSimulateRandomSpots_Result")
 				FolderName          =   os.path.basename(Image_folder)
 				Result_pathway          =   os.path.join("Z:\Result0309\s15_3H1", FolderName, "XTSimulateRandomSpots_Result")
				CreateDirectoryToSaveFiles(Result_pathway)
				AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
				AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
				TotalNumberFile         =           len(AllFilesToBeProcessed)
				logtime('Start bach processing')
				if TotalNumberFile  >   0:
					for vFileName in AllFilesToBeProcessed:
						try:
							vFullFileName = os.path.join(Image_folder, vFileName)
							vImaris.FileOpen(vFullFileName, "")
							logtime('Image processing START image: '+str(vFileName))
							ImageIsEmpty   =           GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName)
							if not ImageIsEmpty:
								FileIndex	               +=  		1
								FileNameList.append(vFileName)
							else:
								print "File "+vFileName+" doesn't have a data set."
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
				ImageIsEmpty            =   True
				if vFileName !=""       :
					vFileName           =   vImaris.GetCurrentFileName()
					vFilePath           =   os.path.dirname(vFileName)
					vFullFileName = os.path.join(vFilePath, vFileName)
					Result_pathway      =   os.path.join(vFilePath, "XTSimulateRandomSpots_Result")
					CreateDirectoryToSaveFiles(Result_pathway)
					vFileName           =   os.path.split(vFileName)[1]
					ImageIsEmpty        =   GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName)
				if not ImageIsEmpty :
					FileIndex	               +=  		1
					FileNameList.append(vFileName)					
				else:
					print "File "+vFileName+" doesn't have a data set."
			logtime('XTSimulateRandomSpots extension done')
			print "All tasks have been completed successfully. \n Resulting files are saved in the folder XTSimulateRandomSpots_Result"
			raw_input("Press Enter to terminate.")
		else:
			tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
			logtime('Extension XTSimulateRandomSpots END')
	except:
		logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
