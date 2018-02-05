# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 13:49:43 2017

@author: Pheonix
"""
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
