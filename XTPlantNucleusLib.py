# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:16:25 2017

@author: Nani
"""

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
        NewColnames                     =       [x+str(y) for x in FileNameList for y in data.columns] 
        DataResult.columns              =       NewColnames
        vPathToSaveTables               =       os.path.join(Result_pathway, FileType)
        DataResult.to_csv(path_or_buf=vPathToSaveTables + ".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
         
def PullSimilarDataIntoOneFile(FileNameList, Result_pathway):  
    AllFiles                =     os.listdir(Result_pathway) #get all files in the Result_pathway directory
    ListTable               =   ["SurfaceIntensity_Ch", "SurfacePositions_Ch", "SurfaceVolume_Ch"]
    AllChannels             =   [DAPIChannel]+SelectedChanelIndex
    for typeTable in ListTable:
        for sp in AllChannels:
            FileType                            =       typeTable+str(sp)
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
#Function to export output in .csv files

def CreateTables(FileNameList, vData, ExpFileName, Result_pathway,ChID, FormatTable, FileIndex):
    if len(vData)>0:
        SaveTable         =       pd.DataFrame(vData)    #Save position table
        if FormatTable   :
            SaveTable     =       SaveTable.T
        if FileNameList != None:
            data                            =       pd.DataFrame(FileNameList)
            SaveTable                      =       pd.concat([SaveTable, data], axis=1)
        vPathToSaveTables = os.path.join(Result_pathway, ExpFileName+str(ChID))
        SaveTable.to_csv(path_or_buf=vPathToSaveTables+str(FileIndex)+".csv", na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    logtime('Data saved- image_'+str(FileIndex))
    
#Funtion to get mean intensity inside nucleus surface for all channels selected by the user
def GetNucleusIntensity(surf, Dyes):
    Dyes                =   [DAPIChannel]   +   SelectedChanelIndex
    vAllStatistics      =     surf.GetStatistics()
    vNames              =     vAllStatistics.mNames
    vValues             =     vAllStatistics.mValues
    IntenityList        =     [float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Sum"]
    IntensitySum       =   [IntenityList[y] for y in Dyes]
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
                mask            =    [maxValue if item>0.7 else item for item in mask]
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

def SegImmunChannel(sel_channel, DiameterLlist, aThreshold,FileIndex, GroupContainer, ContainerSurface, numberIndex):
    vImage                =    vImaris.GetDataSet()
    vSurpassScene         =   vImaris.GetSurpassScene()
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
        vDiameter_XY_index            =    DiameterLlist[diametIndex]
        diametIndex                 +=  1
        vDiameter_Z_index            =    DiameterLlist[diametIndex]
        diametIndex                 +=  1
        #Set the maximum intensity threshold for channel
        aMaximum_intensity_threshold=    (float(aThreshold)/100.0)*Max_intensity[channel_index]
        aRegionsOfInterest          =    None
        aSubtractBackground            =     1
        aSpotFiltersString            =    '"Quality" above automatic threshold'
        aEstimateDiameterXYZ        =     [vDiameter_XY_index, vDiameter_XY_index, vDiameter_Z_index]
        vSpots                        =    vImaris.GetImageProcessing().DetectEllipticSpots(vImage,aRegionsOfInterest,numberIndex-1,aEstimateDiameterXYZ,aSubtractBackground,aSpotFiltersString)
        logtime('Immunostaining signal segemented - channelId: ' + str(channel_index)+' - Image '+str(FileIndex))
        Table                =    vSpots.GetPositionsXYZ()
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
    result                  =   SelectVoxels(MaskSurface, vImage_data, "HighIntensity")
    AddChannel(result, vImage, numberIndex) #CC channel
    vImage.SetChannelName(numberIndex-1,ChannelName)
    Surface                 =   Segment_Surface(vImage, numberIndex-1, SmoothDetail, ObjectName,ObjectDiam)
    GroupContainer.AddChild(Surface, -1)
    vScene.AddChild(GroupContainer, -1)
    NumberOfSurfaces        =    Surface.GetNumberOfSurfaces()
    IntenistyList           =   get_intensity(Dyes, Surface, NumberOfSurfaces)
    positionOfAllItems      =   [Surface.GetCenterOfMass(surfId)[0] for surfId in range(NumberOfSurfaces)]
    SurfVolume              =       GetNucParameters(Surface)[0]
    if type(SurfVolume) is float:
        SurfVolume          =   [SurfVolume]
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
