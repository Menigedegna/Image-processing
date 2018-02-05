# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus, nucleolus and chromocenters into surfaces in DAPI channel,
# 	Create masks using these surfaces,
#   Segments RNA PolII immunostaining signal into spots in Immunostaining channels,
#   Simulate random spots
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
#       <Item name="CountSpotPerShell" icon="Python" tooltip="CountSpotPerShell">
#         <Command>PythonXT::XTCountSpotPerShell(%i)</Command>
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

#==============================================================================
# These functions are required to do spot distribution analysis :
# Create shells for nucleus, for chromocenters and or for nucleolus
# Get masks for each shell
# Use masks to identify voxels inside shells
# Convert spot coordinates into voxels index
# Count number of spots in each mask
#==============================================================================
# Function to generate random spots in nucleus
def simulate3DPointsInSurface(vFactory, groupContainer, aRadius, chID, number_spots, surfMask, VoxelSizeX, VoxelSizeY, VoxelSizeZ):
    """Choose random voxel IDs from mask created"""
    RandomPointVoxelId        =   np.random.choice(range(len(surfMask)), number_spots)
    RandomPointVoxel          =   [surfMask[x] for x in RandomPointVoxelId]
    """Concert voxel ID into X,Y,Z position coordiantes"""
    RandomPosition      =   ConvertCoordianteIntoVoxelId("rev", RandomPointVoxel, VoxelSizeX, VoxelSizeY, VoxelSizeZ)
    SpotName            =   "Random spot - Ch" + str(chID)
    CreateSpots(vFactory, RandomPosition, SpotName, groupContainer, aRadius)
    return RandomPosition, RandomPointVoxel

# Function to create spots
def CreateSpots(vFactory, aPositionsXYZ, SpotName, groupContainer, aRadius):
    vSpot			        =	vFactory.CreateSpots()
    aIndicesT				=	[0.0]*len(aPositionsXYZ)
    aRadii					=	[aRadius]*len(aPositionsXYZ)
    vSpot.Set(aPositionsXYZ,aIndicesT,aRadii)
    vSpot.SetName(SpotName)
    groupContainer.AddChild(vSpot, -1)

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
    label					    =	"Intensity Mean"
    Index			          	=	[i for i, x in enumerate(vNames) if x == label]
    result			            =	[vValues[e] for e in Index]
    IntensityInAllChannel       =   [result[(vChannelI*Number_spot):(vChannelI*Number_spot+Number_spot)] for vChannelI in Dyes]
    NormalisedIntensityList     =   []
    for ch in range(len(MeanIntensitySurfaceList)):
        NormalisedIntensity         =   [round(float(z)/float(MeanIntensitySurfaceList[ch]), 3) for z in IntensityInAllChannel[ch]]
        NormalisedIntensityList.append(NormalisedIntensity)
    return NormalisedIntensityList

# Function to convert spots coordinates into voxel index
def ConvertCoordianteIntoVoxelId(Type, SpotPosition, VoxelSizeX, VoxelSizeY, VoxelSizeZ):
    Result      =   []
    if Type=="Fwd":
        for spot in SpotPosition:
            VoxelID     =   [int(round(x*y, 0)) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
            Result.append(VoxelID)
    else:
        for spot in SpotPosition:
            Coordinates     =   [round(x/y, 3) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
            Result.append(Coordinates)
    return Result

# Function to identify spots inside shells
def GetSpotInsideSurface(radiusList, VolumeShellList, SpotIntensityList, TypeStudy, Coordinatedata, ResultInsideMask, NumberOfShell, vImage, FileIndex, SpotVoxelId,TypeSpot ):
    global ListOfContainers
    vFactory		     =   vImaris.GetFactory()
    vSurpassScene 	     =   vImaris.GetSurpassScene()
    ListIntensityAllSpots   =   []
    ListPointAllChannel     =   []
    NumberOfShell           =   int(NumberOfShell)
    Shellid                 =   range(1,NumberOfShell+1)
    if TypeStudy    !=   "Nucleus":
         ResultInsideMask.reverse()
         VolumeShellList.reverse()
#    VolumeShell             =   [x-y for x,y in zip(VolumeShellList, [0]+VolumeShellList[0:NumberOfShell-1])]
    VolumeShell             =       VolumeShellList
    for sp in range(len(SelectedChanelIndex)):
        ListSpotVoxelId     =   SpotVoxelId[sp]
        PositionTable       =   Coordinatedata[sp]
        TotalNumberSpot     =   len(PositionTable)
        GroupOfShell        =   vFactory.CreateDataContainer()
        GroupOfShell.SetName(TypeSpot +" - Spots' inside "+ TypeStudy+" Ch"+str(SelectedChanelIndex[sp]))
        ListOfContainers.append(GroupOfShell)
        NextShell           =   ResultInsideMask[1:NumberOfShell]
        NextShell.append([None, None, None])
        SpotIntensity       =   SpotIntensityList[sp]
        ListIntensityInAllShells    =   []
        ListPointAllShells  =   []
        for sh in range(NumberOfShell):
            selectedIndex   =   [x for x in range(len(ListSpotVoxelId)) if ListSpotVoxelId[x] in ResultInsideMask[sh] and ListSpotVoxelId[x] not in NextShell[sh]]
            InsidePoints    =   [PositionTable[x] for x in selectedIndex]
            CreateSpots(vFactory, InsidePoints, "Shell - "+str(Shellid[sh]), GroupOfShell, radiusList[sp])
            RelatifNumberSpot   =   round(float(len(InsidePoints))/float(TotalNumberSpot)/float(VolumeShell[sh]), 3)
            ListPointAllShells.extend([RelatifNumberSpot])
            ListSelectedIntensitiesPerChannel =   []
            for IntensityPerChannel in SpotIntensity:
                IntensityList   =   [IntensityPerChannel[y] for y in selectedIndex]
                ListSelectedIntensitiesPerChannel.append(IntensityList)
            ListIntensityInAllShells.append(ListSelectedIntensitiesPerChannel)
        ListIntensityAllSpots.append(ListIntensityInAllShells)
        ListPointAllChannel.append(ListPointAllShells)
        vSurpassScene.AddChild(GroupOfShell, -1)
        logtime("Spots inside shells are selected for ch"+str(SelectedChanelIndex[sp])+" - Image "+str(FileIndex))
    return ListPointAllChannel, ListIntensityAllSpots

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
    NucMask         =   Get_Mask_data(NuucleusSurface) if NuucleusSurface is not None else None
    for SurfaceInstance in ListSurface:
        mask_values                 =   Get_Mask_data(SurfaceInstance)
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
            mask_values   =   MaskX
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

# Function to quantify spots in image
def CountNumberOfSpot(TypeStudy, numberIndex, ListParameters, ImmunoSpotList, ChromocenterSurface, NucleusSurface, NucleolusSurface, NucleusSurfaceVolume,ImmunoNames, FileIndex, Result_pathway, DistanceOptions, DiameterListe):
    global ListOfContainers 
    Dyes                        =   [DAPIChannel]   +   SelectedChanelIndex # Indexes of intensities I save for each spot
    NumberSelectedChannel       =   len(SelectedChanelIndex)
    radiusList                  =   [np.mean(DiameterListe[a:b])/2 for a,b in zip(range(0,NumberSelectedChannel, 2), range(2,NumberSelectedChannel+2, 2))] # radius for spot segmentation for foci channels, these diamters are set by user
    vImage			         =   vImaris.GetDataSet()
    vFactory		         =   vImaris.GetFactory()
    vSurpassScene 		   =   vImaris.GetSurpassScene()
    NumberOfShell               =   int(ListParameters[8])
    MeanIntensitySurfaceList    =   GetNucleusIntensity(NucleusSurface)
    ListMask                    =   [] # Track mask of shells created
    ListOfVolume                =   [] # Track volume of shells created
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
    VoxelSizeX              =   vImageSizeX/ SixeX
    VoxelSizeY              =   vImageSizeY/ SixeY
    VoxelSizeZ              =   vImageSizeZ/ SixeZ
    SurfaceShellID          =   [] 
    logtime('Count number of spots START - image_'+str(FileIndex))
    NucleusVolume       =       GetVolume(NucleusSurface)
    Label                  =    ["Obs", "Sim"] #2 type of spots are created in scene: randome simulated spots and segmented spots of selected channels
    # 		Create shells for the nucleus surface
    # ==============================================================================
    if DistanceOptions[0] and NucleusSurface is not None:
        ListOfShells        =       CreateShells(vFactory, [NucleusSurface], NumberOfShell, "Nucleus", None)
        VolumeList          =       [GetVolume(surf) for surf in ListOfShells]
        SurfVolume          =       [round(float(vol)/float(NucleusVolume),3) for vol in VolumeList]
        ShellMaskList       =       GeteMasks(ListOfShells, "", None)
        SurfaceShellID.extend([0])
        ListMask.append(ShellMaskList)
        ListOfVolume.append(SurfVolume)
        logtime("Shells are created for Nucleus - Image "+str(FileIndex))
    # 		Create shells for the chromocenters and nucleolus if user selected this option
    # ==============================================================================
    SurfaceList             =   [x for x,y in zip([ChromocenterSurface, NucleolusSurface], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
    if len(SurfaceList)>0:
        surfaceAnalised     =   [x for x,y in zip(["Chromocenters", "Nucleolus"], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
        SurfaceName         =   '-'.join(surfaceAnalised)
        ListOfShells        =   CreateShells(vFactory, SurfaceList, NumberOfShell, SurfaceName, ListParameters[10])
        ShellMaskList       =   GeteMasks(ListOfShells, "inverse", NucleusSurface)
        SurfaceShellID.extend([-1])
        ListMask.append(ShellMaskList)
        SurfVolume          =   []
        for mask in ShellMaskList:
            AddChannel(mask, vImage, numberIndex)
            surf            =   Segment_Surface(vImage, numberIndex-1,0.2, "ShellSurface",0) #I create a surface that contains every voxel outside of the chromocenters and inside the nucleus.
            volumeSurf      =   GetVolume(surf)
            SurfVolume.extend([volumeSurf])
        SurfVolume          =       [round(float(X)/float(NucleusVolume),3) for X in SurfVolume]
        ListOfVolume.append(SurfVolume)
        logtime("Shells are created for"+SurfaceName+" - Image "+str(FileIndex))
    # 		Quantify spots in each shell
    # ==============================================================================
    SpotsVoxelList          =   []
    SpotIntensityList       =   []
    NumberOfSpotsList       =   []
    PositionList            =   []
    for sp in range(len(SelectedChanelIndex)):
        PositionTable       =   ImmunoSpotList[sp].GetPositionsXYZ()
        PositionList.append(PositionTable)
        ListSpotVoxelId     =   ConvertCoordianteIntoVoxelId("Fwd", PositionTable, VoxelSizeX, VoxelSizeY, VoxelSizeZ)
        SpotsVoxelList.append(ListSpotVoxelId)
        NumberOfSpots       =   len(PositionTable)
        SpotIntensity       =   get_intensity(Dyes, ImmunoSpotList[sp], NumberOfSpots, MeanIntensitySurfaceList)
        SpotIntensityList.append(SpotIntensity)
        NumberOfSpotsList.extend([NumberOfSpots])
    logtime('Spots coordiantes converted into Voxel IDs_'+str(FileIndex))
    StudyNumberSpotResult               =   []
    StudyIntensityResult                =   []
    StudyNumberSpotResultR              =   []
    StudyIntensityResultR               =   []
    for mk in range(len(ListMask)):
        RandomPointsInsideSurfaceList       =   []
        RandomPointsInsideSurfaceVoxelList  =   []
        ShellMaskList       =   ListMask[mk]
        ResultInsideMask    =   [[]]*NumberOfShell
        for sh in range(len(ShellMaskList)):
            ListSelectedVoxels    =   [[x,y,z] for x in range(vImageSizeX) for y in range(vImageSizeY) for z in range(vImageSizeZ) if ShellMaskList[sh][x][y][z]]
            ResultInsideMask[sh]  =   ResultInsideMask[sh]    +   ListSelectedVoxels
        logtime("Voxels inside shells are selected - Image "+str(FileIndex))
        GroupOfRandomPoints         =	vFactory.CreateDataContainer()
        ListOfContainers.append(GroupOfRandomPoints)
        GroupOfRandomPoints.SetName(" Random spots' inside/outside"+ TypeStudy[mk])
        for sp in range(len(ImmunoSpotList)):
            RandomPointsInsideSurface, RandomPointsInsideSurfaceVoxel   =     simulate3DPointsInSurface(vFactory, GroupOfRandomPoints,  radiusList[sp], SelectedChanelIndex[sp], NumberOfSpotsList[sp], ResultInsideMask[SurfaceShellID[mk]],  VoxelSizeX, VoxelSizeY, VoxelSizeZ)
            RandomPointsInsideSurfaceList.append(RandomPointsInsideSurface)
            RandomPointsInsideSurfaceVoxelList.append(RandomPointsInsideSurfaceVoxel)
            logtime("Simulated random points for channel"+str(sp)+"- Image "+str(FileIndex))
        vSurpassScene.AddChild(GroupOfRandomPoints, -1)
        NumberSpotData, IntensityData        =   GetSpotInsideSurface(radiusList, ListOfVolume[mk], SpotIntensityList, TypeStudy[mk], PositionList, ResultInsideMask, NumberOfShell, vImage, FileIndex, SpotsVoxelList, "Obs")
        NumberSpotDataR, IntensityDataR      =   GetSpotInsideSurface(radiusList, ListOfVolume[mk], SpotIntensityList, TypeStudy[mk], RandomPointsInsideSurfaceList, ResultInsideMask, NumberOfShell, vImage, FileIndex, RandomPointsInsideSurfaceVoxelList, "Sim")
        StudyNumberSpotResult.append(NumberSpotData)
        StudyIntensityResult.append(IntensityData)
        StudyNumberSpotResultR.append(NumberSpotDataR)
        StudyIntensityResultR.append(IntensityDataR)
        logtime('Spots inside each'+TypeStudy[mk]+'shell is quantified'+str(FileIndex))
    # 		Save dataframes containing results
    # ==============================================================================
    ResultNumberOfSpot     =    [StudyNumberSpotResult, StudyNumberSpotResultR]
    ResultIntensity        =    [StudyIntensityResult, StudyIntensityResultR]
    for typeOfSpot in range(len(Label)):
        for mk in range(len(TypeStudy)):
            for sp in range(len(ImmunoSpotList)):
                if len(ResultNumberOfSpot[typeOfSpot][mk][sp])>0:
                    SaveTable         =       pd.DataFrame(ResultNumberOfSpot[typeOfSpot][mk][sp])
                    vPathToSaveTables = os.path.join(Result_pathway, "NumberOfSpot_SP"+str(sp)+"_"+Label[typeOfSpot])
                    SaveTable.to_csv(path_or_buf=vPathToSaveTables +"_ST"+str(mk)+"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
                for sh in range(int(NumberOfShell)):
                    for ch in range(len(Dyes)):
                        if len(ResultIntensity[typeOfSpot][mk][sp][sh][ch])>0:
                            SaveTable         =       pd.DataFrame(ResultIntensity[typeOfSpot][mk][sp][sh][ch])
                            vPathToSaveTables = os.path.join(Result_pathway, "Intensity_SP"+str(sp)+"_SH"+str(sh)+"_CH"+str(Dyes[ch])+"_"+Label[typeOfSpot])
                            SaveTable.to_csv(path_or_buf=vPathToSaveTables +"_ST"+str(mk)+"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    logtime('Data saved- image_'+str(FileIndex))
    logtime('Count number of spots END - image_'+str(FileIndex))

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
def SelectVoxels(mask_values, vImage_data, Type):
    maxValue=getMax(vImage_data)
    Result                  =   []
    for z in range(len(vImage_data)):
        maskY               =   []
        for y in range(len(vImage_data[0])):
            mask_surface	=	mask_values[z][y]
            mask            =	vImage_data[z][y]
            if Type=="HighIntensity":
                mask            =   [round(item/maxValue, 0) for item in mask]
                mask            =	[maxValue if item==1 else item for item in mask]
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

def SegImmunChannel(sel_channel, DiameterLlist, aThreshold,FileIndex, GroupContainer, ContainerSurface, numberIndex):
    vImage			    =	vImaris.GetDataSet()
    vSurpassScene 		=   vImaris.GetSurpassScene()
    ImmunoSpotList      =   []
    SpotNameList        =   []
    diametIndex         =   0
    Max_intensity       =   GetMaxIntensityOfChannels(vSurpassScene)
    for channel_index in sel_channel:
        # Clean channel : select voxels inside nucleus surface only
        vImage_data     =   vImage.GetDataVolumeFloats(channel_index,0)
        MaskSurface     =   Get_Mask_data(ContainerSurface)
        CleanedImage    =   SelectVoxels(MaskSurface, vImage_data, "")
        AddChannel(CleanedImage, vImage, numberIndex)
        vImage.SetChannelName(numberIndex-1,"Immuno_Ch"+str(channel_index))
        vDiameter_XY_index			=	DiameterLlist[diametIndex]
        diametIndex                 +=  1
        vDiameter_Z_index			=	DiameterLlist[diametIndex]
        diametIndex                 +=  1
        #Set the maximum intensity threshold for channel
        aMaximum_intensity_threshold=	(float(aThreshold)/100.0)*Max_intensity[channel_index]
        aRegionsOfInterest          =    None
        aSubtractBackground		    =	 1
        aSpotFiltersString		    =    '"Quality" above automatic threshold'
        aEstimateDiameterXYZ		=	 [vDiameter_XY_index, vDiameter_XY_index, vDiameter_Z_index]
        vSpots				        =	vImaris.GetImageProcessing().DetectEllipticSpots(vImage,aRegionsOfInterest,numberIndex-1,aEstimateDiameterXYZ,aSubtractBackground,aSpotFiltersString)
        logtime('Immunostaining signal segemented - channelId: ' + str(channel_index)+' - Image '+str(FileIndex))
        Table				=	vSpots.GetPositionsXYZ()
        if Table is not None:
            ImmunoSpotList  =   ImmunoSpotList  +   [vSpots]
            SpotName        =   "Immuno_Ch"+str(channel_index)
            vSpots.SetName(SpotName)
            GroupContainer.AddChild(vSpots, -1)
            vSurpassScene.AddChild(GroupContainer, -1)
            SpotNameList    =   SpotNameList  +   [SpotName]
        numberIndex         +=    1
    return ImmunoSpotList, SpotNameList

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

#This function, inside a given surface, selectes voxels containing 50% of the lowest intensity, sets there value to the maximum intensity, and sets the value of the other voxels to 0
def getLowIntensity(mask_values, vImage_data, mask_values2):
    maxValue                =   getMax(vImage_data)
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

#This function first creates a channel by selecting voxel contianing 50% of highest intensity values
def SegHighIntensity(chIndex, surf, ChannelName, ObjectName, numberIndex, GroupContainer, ParametersList):
    vImage                  =   vImaris.GetDataSet()
    vScene                  =   vImaris.GetSurpassScene()
    vImage_data             =   vImage.GetDataVolumeFloats(chIndex,0)
    MaskSurface             =   Get_Mask_data(surf)
    result                  =   SelectVoxels(MaskSurface, vImage_data, "HighIntensity")
    AddChannel(result, vImage, numberIndex) #CC channel
    vImage.SetChannelName(numberIndex-1,ChannelName)
    result              =   Segment_Surface(vImage, numberIndex-1, ParametersList[5], ObjectName,ParametersList[6])#CC segmented
    return result

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
# Function required to concatenate files containing the same information into one table
#==============================================================================
def ConcatenatedTables(FileNameList, Result_pathway,AllFilesWithParam,FileType):
    DataResult                          =       pd.DataFrame()
    for File in AllFilesWithParam:
        SelectedFile                    =       os.path.join(Result_pathway, File)
        data                            =       pd.read_csv(SelectedFile, header='infer',decimal='.')
        DataResult                      =       pd.concat([DataResult, data], axis=1)
        os.remove(SelectedFile)
    if len(DataResult)>0:
        DataResult.columns              =       FileNameList
        vPathToSaveTables               =       os.path.join(Result_pathway, FileType)
        DataResult.to_csv(path_or_buf=vPathToSaveTables + ".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
         
def PullSimilarDataIntoOneFile(FileNameList, Result_pathway, TypeStudy, SelectedChannel, NumberOfShell):  
    Dyes                   =   [DAPIChannel]   +   SelectedChanelIndex # Indexes of intensities I save for each spot
    AllFiles               =     os.listdir(Result_pathway) #get all files in the Result_pathway directory
    Label                  =    ["Obs", "Sim"]
    for mk in range(len(TypeStudy)):
        for typeOfSpot in range(len(Label)):
            for sp in range(len(SelectedChannel)):
                FileType                            =       "NumberOfSpot_SP" + str(sp)+"_"+Label[typeOfSpot]+"_ST"+str(mk)
                AllFilesWithParam                   =       [i for i in AllFiles if FileType in i]
                ConcatenatedTables(FileNameList, Result_pathway,AllFilesWithParam, FileType)
                for sh in range(int(NumberOfShell)):
                    for ch in Dyes:
                        FileType                            =       "Intensity_SP" + str(sp)+"_SH"+str(sh)+"_CH"+str(ch)+"_"+Label[typeOfSpot]+"_ST"+str(mk)
                        AllFilesWithParam                   =       [i for i in AllFiles if FileType in i]
                        ConcatenatedTables(FileNameList, Result_pathway,AllFilesWithParam, FileType)
                        
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
def GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing):
    global vImaris
    global numberIndex
    global SelectedChanelIndex, DAPIChannel
    global UserSetDiameters
    global ListOfContainers
    ListOfContainers    =   []   # To keep track of all the containers I will create in scene, so that I can remove them after saving image
    vImage              =   	vImaris.GetDataSet()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        if FileIndex    == 1:
            UserSetDiameters=[]
            SelectedChanelIndex=[]
            SelectedChanelIndex, DAPIChannel= Ask_user(numberIndex)
            if len(SelectedChanelIndex)>0:
                AskUserDiameters(SelectedChanelIndex)
            if len(SelectedChanelIndex)==0 and len(UserSetDiameters)==0:
                return 1,1
        numberIndex     +=1
        ChromocenterSurface =   NucleusSurface  =   NucleolusSurface    =   None
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        vImage                     =    vImaris.GetDataSet()
        if DoSegmentation:
            logtime('Segmentation START - image_'+str(FileIndex))
            vScene                     =    vImaris.GetSurpassScene()
            vFactory		           =	vImaris.GetFactory().CreateFactory()
            GroupOfObjects	           =	vFactory.CreateDataContainer()
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
            SurfP=Segment_Surface(vImage, DAPIChannel, ParametersList[4], "Rougher Nucleus segmentation",0)
            NucleusSurface=SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
            logtime('Nucleus surface segmented - Image '+str(FileIndex))
            if DistanceOptions[2]:
                ChromocenterSurface          =  SegHighIntensity(DAPIChannel, SurfP,"High DAPI intensity",
                                                                 "Chromocenters",
                                                                 numberIndex, GroupOfObjects, ParametersList)#Chromocenters are segmented
                GroupOfObjects.AddChild(ChromocenterSurface, -1)
                vScene.AddChild(GroupOfObjects, -1)
                logtime('Chromocenter surfaces segmented - Image '+str(FileIndex))
                numberIndex     +=  1
            ImmunoSignalList     =    []
            ImmunoNames              =    []
            ImmunoSignalList, ImmunoNames     =    SegImmunChannel(SelectedChanelIndex, UserSetDiameters, ParametersList[7],FileIndex, GroupOfObjects, NucleusSurface, numberIndex) #Immunostaining channels are segmented
            logtime('Segmentation END - image_'+str(FileIndex))
            numberIndex =   numberIndex+len(SelectedChanelIndex)
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface,ImmunoNames,ImmunoSignalList  = GetSegmentedObjects(DistanceOptions, FileIndex)
        NucleusSurfaceVolume    =   GetVolume(NucleusSurface)
        CountNumberOfSpot(TypeStudy, numberIndex, ParametersList, ImmunoSignalList, ChromocenterSurface, NucleusSurface, NucleolusSurface, NucleusSurfaceVolume,ImmunoNames, FileIndex, Result_pathway, DistanceOptions, UserSetDiameters)
        if DoSegmentation:
            vPath = os.path.join(Result_pathway, vFileName+".ims")
            vImaris.FileSave(vPath, "")
    if len(ListOfContainers)>0:
        RemoveObjectsCreated(vScene, ListOfContainers)
    return vImage is None,ImmunoNames

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
    Messge="Please select the DAPI channel: \n Please only choose one channel."
    PopUpMessage(OPTIONS, Messge)
    DAPIChannel	=	[i for i, x in enumerate(User_selection) if x == 1][0]
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
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTCountSpotPerShell_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTCountSpotPerShell(aImarisId):
	logging.basicConfig(level=logging.DEBUG, filename= "log[XTCountSpotPerShell].log")
	try:
		#Declare global variables
		global gLasttime
		global vImaris
		global SelectedChanelIndex
		gLasttime               =    None
		logtime('Extension XTCountSpotPerShell START')
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
				Result_pathway          =           os.path.join(Image_folder, "XTCountSpotPerShell_Result")
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
						ImageIsEmpty,ImmunoNames    =           GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing)
						if ImageIsEmpty !=  ImmunoNames and not ImageIsEmpty :
							FileIndex	               +=  		1
							FileNameList.append(vFileName)
						if ImageIsEmpty ==  ImmunoNames:
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
					vFilePath           =   os.path.dirname(vFileName)
					Result_pathway      =   os.path.join(vFilePath, "XTCountSpotPerShell_Result")
					CreateDirectoryToSaveFiles(Result_pathway)
					vFileName           =   os.path.split(vFileName)[1]
					ImageIsEmpty,ImmunoNames        =   GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing)
				if ImageIsEmpty !=  ImmunoNames and not ImageIsEmpty :
					FileIndex	               +=  		1
					FileNameList.append(vFileName)					
				if ImageIsEmpty !=  ImmunoNames and vFileName =="" or ImageIsEmpty :
					tkMessageBox.showinfo(title="Alert", message="No image is detected. \n Please open an image and select on 'CountSpotPerShell' again.")
					quit()
				if ImageIsEmpty ==  ImmunoNames:
					tkMessageBox.showinfo(title="Alert", message="User needs to select channels and set diameter values for spot segmentation.")
					quit()
			if FileIndex>1:
				logtime('Files saved and organised')
				PullSimilarDataIntoOneFile(FileNameList, Result_pathway, TypeStudy, SelectedChanelIndex, ParametersList[8])
			logtime('XTCountSpotPerShell extension done')
			print "All tasks have been completed successfully. \n Resulting files are saved in the folder XTCountSpotPerShell_Result"
			raw_input("Press Enter to terminate.")
		else:
			tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
			logtime('Extension XTCountSpotPerShell END')
	except:
		logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
