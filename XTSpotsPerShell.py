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
# Creator: Mariamawit S. Ashenafi, UZH
# Published on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTSpotsPerShell" icon="Python" tooltip="XTSpotsPerShell">
#         <Command>PythonXT::XTSpotsPerShell(%i)</Command>
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

#==============================================================================
# These functions are required to do spot distribution analysis :
# Create shells for nucleus, for chromocenters and or for nucleolus
# Get masks for each shell
# Use masks to identify voxels inside shells
# Convert spot coordinates into voxels index
# Count number of spots in each mask
#==============================================================================

def simulate3DPointsInSurface(vFactory, groupContainer, aRadius, chID, number_spots, surfMask, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMaxX, vExtentMaxY, vExtentMaxZ):
    """Choose random voxel IDs from mask created"""
    RandomPointVoxelId        =   np.random.choice(range(len(surfMask)), number_spots)
    RandomPointVoxel          =   [surfMask[x] for x in RandomPointVoxelId]
    """Concert voxel ID into X,Y,Z position coordiantes"""
    RandomPosition      =   ConvertCoordianteIntoVoxelId("Rev", RandomPointVoxel, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMaxX, vExtentMaxY, vExtentMaxZ)
    SpotName            =   "Random spot - Ch" + str(chID)
    vSpot               =   CreateSpots(vFactory, RandomPosition, SpotName, groupContainer, aRadius)
    return RandomPosition

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
def ConvertCoordianteIntoVoxelId(Type, SpotPosition, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMinX, vExtentMinY, vExtentMinZ):
    Result      =   []
    if Type=="Fwd":
        for spot in SpotPosition:
            VoxelID     =   [int(round(x*y, 0)) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
            Result.append(VoxelID)
    else:
        for spot in SpotPosition:
            Coordinates     =   [round(x*y, 3) for x,y in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ])]
#            Coordinates     =   [round(x*y+z, 3) for x,y,z in zip(spot, [VoxelSizeX, VoxelSizeY, VoxelSizeZ],[vExtentMinX, vExtentMinY, vExtentMinZ])]
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

def GeteMasks(ListSurface, Type, NuucleusSurface):
    ListSurfaceMasks   =   []
    NucMask         =   SN.Get_Mask_data(NuucleusSurface) if NuucleusSurface is not None else None
    for countShell in range(len(ListSurface)):
        mask_values                 =   SN.Get_Mask_data(ListSurface[countShell])
        MaskPrevShell               =   SN.Get_Mask_data(ListSurface[countShell-1]) if countShell>0 else None
        if Type=="inverse" or countShell>0:
            MaskX       =   []
            for x in range(len(mask_values)):
                MaskY       =   []
                for y in range(len(mask_values[x])):
                    MaskZ   =   mask_values[x][y]
                    MaskNuc =   NucMask[x][y]
                    if Type=="inverse" and countShell==0:
                        MaskZ   =   [b if a==0 else 0 for a,b in zip(MaskZ, MaskNuc)]
                    if countShell>0:
                        MaskPrev     = MaskPrevShell[x][y]
                        MaskPrev   =   [1 if a==0 else 0 for a in MaskPrev]
                        MaskZ   =   [b if a==0 else 0 for a,b in zip(MaskZ, MaskPrev)]                    
                    MaskY.append(MaskZ)
                MaskX.append(MaskY)
            mask_values=MaskX
        ListSurfaceMasks.append(mask_values)
    return ListSurfaceMasks

def transformVertices(x, delta, centerOfMass):
    result=x+delta if (x-centerOfMass)>=0 else x-delta
    return result

def CreateShells(vFactory, SurfaceList, number_of_shell, SurfaceName, MaxFactor):
    global ListOfContainers
    vSurpassScene 		=   vImaris.GetSurpassScene()
    ListOfShells        =   []
    GroupOfShell        =	vFactory.CreateDataContainer()
    GroupOfShell.SetName(SurfaceName + ' shells')
    ListOfContainers.append(GroupOfShell)
    vTimeIndex			=	0
    number_of_shell     =   int(number_of_shell)
    delta=0.2 #the difference between shells is 100nm
#    FactorIncrease      =   0
#    if SurfaceName  !=  "Nucleus":
#        Stepper           =   (float(MaxFactor)-1.0)/float(number_of_shell) # increasing factor
    for shell in range(number_of_shell, -1, -1): # for each shell to be created
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
                if SurfaceName  ==  "Nucleus" and number_of_shell>0:
                    vNewVertices        =   [[(x[0]-vCenterOfMass[0])*(float(shell)/number_of_shell)+vCenterOfMass[0], (x[1]-vCenterOfMass[1])*(float(shell)/number_of_shell)+vCenterOfMass[1], (x[2]-vCenterOfMass[2])*(float(shell)/number_of_shell)+vCenterOfMass[2]] for x in vVertices]
                else:
#                    FactorIncrease      +=  Stepper
#                    vNewVertices        =   [[(x[0]-vCenterOfMass[0])*((number_of_shell+float(shell))/number_of_shell)+vCenterOfMass[0],(x[1]-vCenterOfMass[1])*((number_of_shell+float(shell))/number_of_shell)+vCenterOfMass[1], (x[2]-vCenterOfMass[2])*((number_of_shell+float(shell))/number_of_shell)+vCenterOfMass[2]] for x in vVertices]
                    vNewVertices        =   [[transformVertices(coord, shell*delta, centerOfMass) for coord, centerOfMass in zip(x,vCenterOfMass)] for x in vVertices]
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
    radiusList                  =   [np.mean(DiameterListe[a:b])/2 for a,b in zip(range(0,len(DiameterListe), 2), range(2,len(DiameterListe)+2, 2))] # radius for spot segmentation for foci channels, these diamters are set by user
    vImage			         =   vImaris.GetDataSet()
    vFactory		         =   vImaris.GetFactory()
    vSurpassScene 		   =   vImaris.GetSurpassScene()
    NumberOfShell               =   3
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
    logtime('Count number of spots START - image_'+str(FileIndex))
    Label                  =    ["Obs", "Sim"] #2 type of spots are created in scene: randome simulated spots and segmented spots of selected channels
    # 		Get mask of nucleus surface
    # ==============================================================================
    if DistanceOptions[0] and NucleusSurface is not None:
        ShellMaskList       =       GeteMasks([NucleusSurface], "", None)
        ListMask.extend(ShellMaskList)
        logtime("Shells are created for Nucleus - Image "+str(FileIndex))
    # 		Create mask excluding chromocenters and nucleolus, mask2, if user selected this option
    # ==============================================================================
    SurfaceList             =   [x for x,y in zip([ChromocenterSurface, NucleolusSurface], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
    if len(SurfaceList)>0:
        surfaceAnalised     =   [x for x,y in zip(["Chromocenters", "Nucleolus"], [DistanceOptions[2], DistanceOptions[1]]) if y ==1]
        SurfaceName         =   '-'.join(surfaceAnalised)
        ListOfShells        =   CreateShells(vFactory, SurfaceList, NumberOfShell, SurfaceName, ListParameters[10])
        ShellMaskList       =   GeteMasks(ListOfShells, "inverse", NucleusSurface)
        ListMask.extend(ShellMaskList)
        logtime("Shells are created for"+SurfaceName+" - Image "+str(FileIndex))
    # 		Quantify spots and intensities
    # ==============================================================================
    StudyNumberSpotResult               =   [] # get number of observed spots

    GroupOfRandomPoints                 =   vFactory.CreateDataContainer() # Container for simulated spots
    GroupOfRandomPoints.SetName("Random spots")
    ListOfContainers.append(GroupOfRandomPoints)
    GroupOfShellSpots                 =   vFactory.CreateDataContainer() # Container for simulated spots
    GroupOfShellSpots.SetName("Shell observed spots")
    ListOfContainers.append(GroupOfShellSpots)
    NumberOfSimulation      =   1
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for sp in range(len(SelectedChanelIndex)):
            PositionTable       =   ImmunoSpotList[sp].GetPositionsXYZ()
            if  PositionTable is not None:
#==============================================================================
#                 Export spot positions
#==============================================================================
                Positiondf         =       pd.DataFrame(PositionTable)
                vPathToSaveTables = os.path.join(Result_pathway, "Position_SP"+str(sp)+"_Obs")
                Positiondf.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.') 
#==============================================================================
#                 Get intensity inside each observed spots
#==============================================================================
                NumberOfSpots       =   len(PositionTable)
                SpotIntensity       =   get_intensity(Dyes, ImmunoSpotList[sp], NumberOfSpots, MeanIntensitySurfaceList)
        #        StudyIntensityResult.append(SpotIntensity)
                StudyNumberSpotResult.append(NumberOfSpots)
        #        Export data
                SaveTable         =       pd.DataFrame(SpotIntensity).T
                vPathToSaveTables = os.path.join(Result_pathway, "Intensity_SP"+str(sp)+"_Obs")
                SaveTable.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')  
                logtime('Got number of observed spots and intensities for ChID'+str(SelectedChanelIndex[sp])+"-"+str(FileIndex))
        #   Create random spots in nucleus mask and mask2
        #        MaskList=[]
                for mk in range(len(ListMask)):
                    ShellMaskList           =   ListMask[mk]
                    ListSelectedVoxels    =   [[x,y,z] for x in range(vImageSizeX) for y in range(vImageSizeY) for z in range(vImageSizeZ) if ShellMaskList[x][y][z]]                   
#                    CoordInShell      =   ConvertCoordianteIntoVoxelId("Rev", ListSelectedVoxels, VoxelSizeX, VoxelSizeY, VoxelSizeZ,  vExtentMaxX, vExtentMaxY, vExtentMaxZ)
#                    SelectedSpots       =   [x for x in PositionTable if x in  CoordInShell]
#                    vSpot               =   CreateSpots(vFactory, SelectedSpots, "Spot"+str(sp)+"Shell"+str(mk), GroupOfShellSpots, radiusList[sp])
                    SN.AddChannel(ShellMaskList, vImage, numberIndex)
                    ImmunoSignalList, ImmunoNames     =   SegImmunChannel([numberIndex-1], UserSetDiameters, ListParameters[7],FileIndex, GroupOfShellSpots, NucleusSurface, numberIndex)
                    vSurpassScene.AddChild(GroupOfShellSpots, -1)
#==============================================================================
#         Simulate random spots
#==============================================================================
                    with concurrent.futures.ProcessPoolExecutor() as executor:
                        for sim in range(NumberOfSimulation):
                            RandomSpot          =   simulate3DPointsInSurface(vFactory, GroupOfRandomPoints,  radiusList[sp], SelectedChanelIndex[sp], NumberOfSpots, ListSelectedVoxels,  VoxelSizeX, VoxelSizeY, VoxelSizeZ, vExtentMaxX, vExtentMaxY, vExtentMaxZ)
                            PositionTable         =       pd.DataFrame(RandomSpot)
                            vPathToSaveTables = os.path.join(Result_pathway, "Position_SP"+str(sp)+"_ST"+str(mk)+"_Sml"+str(sim)+"_Sim")
                            PositionTable.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.') 
                logtime('Got intensities simulated spots for ChID'+str(SelectedChanelIndex[sp])+"-"+str(FileIndex))
    vSurpassScene.AddChild(GroupOfRandomPoints, -1)
    # 		Save dataframes containing results
    # ==============================================================================
    SaveTable         =       pd.DataFrame(StudyNumberSpotResult)
    vPathToSaveTables = os.path.join(Result_pathway, "NumberOfSpot_SP_Obs_")
    SaveTable.to_csv(path_or_buf=vPathToSaveTables +"_"+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.') 
    logtime('Data saved- image_'+str(FileIndex))
    logtime('Count number of spots END - image_'+str(FileIndex))
                      
#==============================================================================
# This function: removes all objects created in scene
#==============================================================================
def RemoveObjectsCreated(vScene, ListOfContainers):
    for i in ListOfContainers:
        vScene.RemoveChild(i)
        
def SegImmunChannel(sel_channel, DiameterLlist, aThreshold,FileIndex, GroupContainer, ContainerSurface, numberIndex):
    vImage			    =	vImaris.GetDataSet()
    vSurpassScene 		=   vImaris.GetSurpassScene()
    ImmunoSpotList      =   []
    SpotNameList        =   []
    diametIndex         =   0
    Max_intensity       =   SN.GetMaxIntensityOfChannels(vSurpassScene)
    for channel_index in sel_channel:
        # Clean channel : select voxels inside nucleus surface only
        vImage_data     =   vImage.GetDataVolumeFloats(channel_index,0)
        MaskSurface     =   SN.Get_Mask_data(ContainerSurface)
        CleanedImage    =   SN.SelectVoxels(MaskSurface, vImage_data, "", ContainerSurface, channel_index)
        SN.AddChannel(CleanedImage, vImage, numberIndex)
        vImage.SetChannelName(numberIndex-1,"Immuno_Ch"+str(channel_index))
        ChannelColor   =   vImage.GetChannelColorRGBA(channel_index)
        vImage.SetChannelColorRGBA (numberIndex-1, ChannelColor)
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
# This function:
# Segments chanenls into surface (nucleus, nucleolus and chromocenters) or spots (RNA PolII foci),
# Create masks with the surfaces created
# Count the number of spots in each mask
# And saves results
#==============================================================================
def GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName):
    global vImaris
    global numberIndex
    global SelectedChanelIndex, DAPIChannel
    global UserSetDiameters
    global ListOfContainers
    ListOfContainers    =   []   # To keep track of all the containers I will create in scene, so that I can remove them after saving image
    ImmunoSignalList     =    []
    ImmunoNames              =    []
    ChromocenterSurface =   NucleusSurface  =   NucleolusSurface    =   None
    vImage              =   	vImaris.GetDataSet()
    vScene                     =    vImaris.GetSurpassScene()
    vFactory		           =	vImaris.GetFactory().CreateFactory()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        if FileIndex    == 1:
            UserSetDiameters=[]
            SelectedChanelIndex=[]
            SelectedChanelIndex= Ask_user(numberIndex)
            if len(SelectedChanelIndex)>0:
                AskUserDiameters(SelectedChanelIndex)
            if len(SelectedChanelIndex)==0 and len(UserSetDiameters)==0:
                return 1,1
        numberIndex     +=1
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        GroupOfObjects	           =	vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
        ListOfContainers.append(GroupOfObjects)
#==============================================================================
# SEGMENT SURFACES
#==============================================================================
        vImage                     =    vImaris.GetDataSet()
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
            IsImageCorrect, NucleusSurface,ChromocenterSurface,  NucleolusSurface, DAPIChannel=SN.SegmentAndGetFeatures(vImage, FileIndex, Result_pathway, vFileName, DistanceOptions, ParametersList, BatchProcessing, vFullFileName, True)
#==============================================================================
#  SEGMENT FOCI
#==============================================================================
            if NucleusSurface.GetNumberOfSurfaces()>0:
                numberIndex =   numberIndex+2
                ImmunoSignalList, ImmunoNames     =   SegImmunChannel(SelectedChanelIndex, UserSetDiameters, ParametersList[7],FileIndex, GroupOfObjects, NucleusSurface, numberIndex)
                numberIndex =   numberIndex+len(SelectedChanelIndex)
            logtime('Segmentation END - image_'+str(FileIndex))
        else:
            NucleusSurface, NucleolusSurface, ChromocenterSurface,ImmunoNames,ImmunoSignalList  = GetSegmentedObjects(DistanceOptions, FileIndex)
#==============================================================================
# GET FOCI FEATURES 
#==============================================================================
        if NucleusSurface.GetNumberOfSurfaces()>0 and ImmunoSignalList     !=    []:
            NucleusSurfaceVolume         =       SN.GetStat(NucleusSurface, "Volume")
            CountNumberOfSpot(TypeStudy, numberIndex, ParametersList, ImmunoSignalList, ChromocenterSurface, NucleusSurface, NucleolusSurface, NucleusSurfaceVolume,ImmunoNames, FileIndex, Result_pathway, DistanceOptions, UserSetDiameters)
#            vPath = os.path.join(Result_pathway, vFileName+".ims")
#            vImaris.FileSave(vPath, "")
#            NucleusSurface.SetVisible(0) #to save snaphot of the the nucleolus and the chromocenters           
#            vPathToSaveTables = os.path.join(Result_pathway, "Snapshot_"+vFileName+".tif")
#            vImaris.SaveSnapShot(vPathToSaveTables)
    else:
        print ("No image detected in file: "+vFileName)
        quit()
    if len(ListOfContainers)>0 and BatchProcessing:
        RemoveObjectsCreated(vScene, ListOfContainers)
#    os.remove(vFullFileName)
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
    Messge                          =       "Please select the RNA PolII immunostaining channels: \n You can choose multiple channels."
    PopUpMessage(OPTIONS, Messge)
    SelectedChanelIndex	            =	    [i for i, x in enumerate(User_selection) if x == 1]
    return SelectedChanelIndex

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
            print ("You must enter a number in all the fields")
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
        tkMessageBox.showinfo(title="Alert", message="Please save the folder 'XTSpotsPerShell_Result' under another name first!")
        quit()
    else:
        os.makedirs(Result_pathway)

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTSpotsPerShell(aImarisId):
	logging.basicConfig(level=logging.DEBUG, filename= "log[XTSpotsPerShell].log")
	try:
		#Declare global variables
		global gLasttime
		global vImaris
		global SelectedChanelIndex
		gLasttime               =    None
		logtime('Extension XTSpotsPerShell START')
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
				FolderName          =   os.path.basename(Image_folder)
				Result_pathway          =   os.path.join("Z:\Result0309\s16_3H1", FolderName, "XTSpotsPerShell_Result")
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
        						ImageIsEmpty,ImmunoNames    =           GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName)
        						if ImageIsEmpty !=  ImmunoNames and not ImageIsEmpty :
        							FileIndex	               +=  		1
        							FileNameList.append(vFileName)
        						if ImageIsEmpty ==  ImmunoNames:
        							tkMessageBox.showinfo(title="Alert", message="User needs to select channels and set diameter values for spot segmentation.")
        							quit()
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
					Result_pathway      =   os.path.join(vFilePath, "XTSpotsPerShell_Result")
					CreateDirectoryToSaveFiles(Result_pathway)
					vFileName           =   os.path.split(vFileName)[1]
					ImageIsEmpty,ImmunoNames        =   GetImageFeatures(TypeStudy, FileIndex, Result_pathway, vFileName, DoSegmentation, DistanceOptions, ParametersList, BatchProcessing, vFullFileName)
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
			logtime('XTCountSpotPerShell extension done')
			print ("All tasks have been completed successfully. \n Resulting files are saved in the folder XTCountSpotPerShell_Result")
			raw_input("Press Enter to terminate.")
		else:
			tkMessageBox.showinfo(title="Alert", message="Imaris application is not found!")
			logtime('Extension XTSpotsPerShell END')
	except:
		logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
