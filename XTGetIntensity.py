# -*- coding: utf-8 -*-
#
#
#==============================================================================
# Objectives of this PythonXT for Imaris:
#   Segments nucleus into a surface in DAPI channel,
#   Get mean and sum DAPI intensity inside surface
#   Get volume of surface
#   Cluster nucleus into ploidy
#   Exports result into .csv tables
# Note: This script is calibrated for 3D images of plant nuclei conterstained in DAPI obtained using Leica TCS SP8
# Creator: Mariamawit S. Ashenafi, Célia Baroux, UZH
# Published on 23.01.2016
#==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTGetIntensity" icon="Python" tooltip="GetIntensity">
#         <Command>PythonXT::XTGetIntensity(%i)</Command>
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
#==============================================================================
# The following python packages are only installed on my home directory on the virtual machines
# TODO: Tell user to add any additional directoy containing packages in Preferences/Custom Tool
# Note this adds more time to load  all libs required
# Keep this option for zmb users or edit it out before publishing 
#==============================================================================
import sys #TODO: edit out?
sys.path.insert(1,'H:\Python') #TODO: edit out?
import scipy #TODO: edit out?
from sklearn.cluster import KMeans #both numpy and scipy needs to be imported before importing sklearns 
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
    Messge="Please select the DAPI channel: \n Please only choose one channel."
    PopUpMessage(OPTIONS, Messge)
    DAPIChannel	=	[i for i, x in enumerate(User_selection) if x == 1][0]
    return DAPIChannel

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
def GetNucleusIntensity(surf):
    vAllStatistics 		= 	surf.GetStatistics()
    vNames	       		= 	vAllStatistics.mNames
    vValues        		= 	vAllStatistics.mValues
    IntenityList        = 	[float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Mean"]
    IntensityMean       =   IntenityList[DAPIChannel]
    IntenityList        =   [float(vValues[a]) for a, x in enumerate(vNames) if x == "Intensity Sum"]
    IntensitySum       =   IntenityList[DAPIChannel]
    return IntensityMean, IntensitySum
    
    
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
def fetch(entries, root):
    global SmothingFactor
    text = entries.get()
    try:
        x = float(text)
    except ValueError:
        print "You must enter a number in field"
        return
    root.destroy()
    SmothingFactor = x


def makeform(root):
    field = "nucleus segmentation"
    row = Frame(root)
    lab = Label(row, width=30, text=field, anchor='w')
    ent = Entry(row)
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ent.pack(side=RIGHT, expand=YES, fill=X)
    return ent


def AskUserSmoothingFactor():
    root = Tk()
    label_text = "Set smooth surface detail for "
    row = Frame(root)
    lab = Label(row, width=30, text=label_text, anchor='w')
    row.pack(side=TOP, fill=X, padx=5, pady=5)
    lab.pack(side=LEFT)
    ents = makeform(root)
    root.bind('<Return>', (lambda event, e=ents: fetch(e, root)))
    b1 = Button(root, fg="darkgreen", text='Submit',
                command=(lambda e=ents: fetch(e, root)))
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
        
#==============================================================================
# This function:
# Segments DAPI chanenl into nucleus surface and gets the mean and max DAPI intensity inside nucleus and the volume of the nucleus
# All data is then registered into IntensityMax.csv,IntensityMean.csv and VolumeNuclei.csv
#==============================================================================
def GetImageFeatures(FileIndex, Result_pathway, vFileName, vFullFileName):
    global vImaris
    global DAPIChannel
    global SmothingFactor
    vImage               =   	vImaris.GetDataSet()
    if vImage is not None:
        numberIndex         = vImage.GetSizeC()
        if FileIndex    == 1:
            DAPIChannel= Ask_user(numberIndex)
            AskUserSmoothingFactor()
        NucleusSurface     =   None
        date			=	str(datetime.datetime.now()).split(" ")[0]
        date			+=	" 00:00:00"
        vImage.SetTimePoint(0, date) # time point needs to be set for each image to avoid creating different time frames when closing and opening a new image.
        vImage                    =    vImaris.GetDataSet()
        logtime('Segmentation START - image_'+str(FileIndex))
        vScene                    =    vImaris.GetSurpassScene()
        vFactory		           =	vImaris.GetFactory().CreateFactory()
        GroupOfObjects	           =	vFactory.CreateDataContainer()
        GroupOfObjects.SetName('Segmented objects')
        SurfP=Segment_Surface(vImage, DAPIChannel, SmothingFactor, "Rougher Nucleus segmentation",0)
        NucleusSurface=SelectSurface(vScene, SurfP, "Nucleus", GroupOfObjects) #Nucleus is segmented
        logtime('Nucleus surface segmented - Image '+str(FileIndex))
        if NucleusSurface is not None: 
            NucMask             =       Get_Mask_data(NucleusSurface) 
            NucMask             =       pd.DataFrame(NucMask)
            vPathToSaveTables = os.path.join(Result_pathway, "Mask_"+str(FileIndex)+".csv")
            NucMask.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
        vScene.RemoveChild(GroupOfObjects) #remove containter from scene
        os.remove(vFullFileName)
    return vImage is None


def BarplotFigure(ClusterIndex, SumIntensityList,Result_pathway):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ## the data
    UniqueValuesClusters        = list(set(ClusterIndex))
    NumberOfImages              = len(SumIntensityList)
    NumberOfNuclei              = [float(ClusterIndex.count(n)) for n in UniqueValuesClusters]
    PercentageOfNucleiPerCluster    = [(float(ClusterIndex.count(n)) / float(NumberOfImages)) * 100.0 for n in
                                UniqueValuesClusters]
    AverageIntensity            = []
    IntensityStd                = []
    for clusterID in UniqueValuesClusters:
        SelectedValuesIntensity = [SumIntensityList[x] for x in range(NumberOfImages) if ClusterIndex[x] == clusterID]
        AverageIntensity.append(np.mean(SelectedValuesIntensity))
        IntensityStd.append(np.std(SelectedValuesIntensity))
    #Order clusters according to DAPI intensity average value 
    SortedClusterIndexValue            =   [x for (y,x) in sorted(zip(AverageIntensity,UniqueValuesClusters))]
    n=0
    for i in SortedClusterIndexValue:
        UniqueValuesClusters[i]                    =   "P"+str(n)
        n+=1
    ind = np.array(range(len(UniqueValuesClusters)))  # the x locations for the groups
    width = 0.35  # the width of the bars
    ## the bars
    rects1 = ax2.bar(ind, AverageIntensity, width,
                     color='olive',
                     yerr=IntensityStd,
                     error_kw=dict(elinewidth=2, ecolor='olive'))
    rects2 = ax.bar(ind + width, PercentageOfNucleiPerCluster, width,
                    color='mediumslateblue')
    # axes and labels
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0.0, 100.0)
    ax.set_ylabel('Percentage of nuclei', color="mediumslateblue")
    ax.tick_params(axis='y', colors='mediumslateblue')
    ax.tick_params(axis='x', colors='white')
    ax2.set_ylabel('Average DAPI intensity sum', color="olive")
    ax2.tick_params(axis='y', colors='olive')
    ax.set_title('DAPI intensity sum and percentage of nuclei by cluster', color="white")
#    xTickMarks = ['Cluster' + str(i) for i in range(1, len(UniqueValuesClusters) + 1)]
    xTickMarks = UniqueValuesClusters
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ## add a legend
    ax.legend((rects1[0], rects2[0]), ('Intensity sum average', 'Number of nuclei'), fontsize=8)
    PlotFile = os.path.join(Result_pathway, "PloidyPlot.png") 
    fig.savefig(PlotFile, transparent=True)
    plt.close(fig)
    dataR       =       pd.DataFrame({"Cluster": UniqueValuesClusters, "SumIntensity Average":  AverageIntensity, 'SumIntensity STD': IntensityStd, 'Number of Nuclei': NumberOfNuclei})
    vPathToSaveTables = os.path.join(Result_pathway,"ClusterIntensityAverage.csv")
    dataR.to_csv(path_or_buf=vPathToSaveTables, na_rep='', float_format=None, columns=None, header=True, index=False, decimal='.')
    return UniqueValuesClusters

#==============================================================================
# Main function:
# Connects to Imaris and get image
# Process images
#==============================================================================
def XTGetIntensity(aImarisId):
    logging.basicConfig(level=logging.DEBUG, filename= "log[XTGetIntensity].log")
    try:
        #Declare global variables
        global gLasttime
        global vImaris
        gLasttime               =    None
        logtime('Extension XTGetIntensity START')
        FileNameList            =   []
        VolumeList               =   []
        MeanIntensityList               =   []
        SumIntensityList               =   []
        #            Step1: Connect to Imaris
        #==============================================================================
        vImarisLib			=	ImarisLib.ImarisLib()
        # Get an imaris object with id aImarisId
        vImaris             =   vImarisLib.GetApplication(aImarisId)
        # """BEGIN LOOP OVER ALL IMAGES"""
        # Open File and get filename
        if vImaris is not None :
            logtime('Connected to Imaris')
            FileIndex                       =       1 #This  variable is used to count the number of files analysed
            # 		Step2: Here the user is asked to set the path to the folder containing image to be analysed
            #==============================================================================
            root1				 =	        Tk()
            Image_folder	       =	        tkFileDialog.askdirectory(parent=root1, initialdir="/",title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
            root1.destroy()
            Result_pathway          =           os.path.join(Image_folder, "XTGetIntensity_Result")
            CreateDirectoryToSaveFiles(Result_pathway)
            AllFilesInDirectory     =           os.listdir(Image_folder) #get all files in the Image_folder directory
            logtime('Get all files')
            AllFilesToBeProcessed   =           [i for i in AllFilesInDirectory if i.endswith('.ims') or i.endswith('.ics')] #select files with .ims or .ics extensions
            TotalNumberFile         =           len(AllFilesToBeProcessed)
            logtime('Start bach processing')
            NumberFISHList                   =       []
            if TotalNumberFile  >   0:
                for vFileName in AllFilesToBeProcessed:
                    vFullFileName   =           os.path.join(Image_folder, vFileName)
                    vImaris.FileOpen(vFullFileName, "")
                    ImageIsEmpty   =           GetImageFeatures(FileIndex, Result_pathway, vFileName, vFullFileName)
                    if not ImageIsEmpty :
                        FileIndex	               +=  		1
                        FileNameList.append(vFileName)
                #Ploidy analysis will be based on Sum of DAPI intensity and nucleus volume fit
            else:
                tkMessageBox.showinfo(title="Alert", message="There is no .ims or .ics file detected in the selected folder.")
                quit()
            logtime('Extension XTGetIntensity END')
            print "All tasks have been completed successfully. \n Resulting files are saved in the folder XTGetIntensity_Result"
            raw_input("Press Enter to terminate.")
        else:
            messagebox.showinfo(title="Alert", message="Imaris application is not found!")
            logtime('Extension XTGetIntensity END')
    except:
        logging.exception("Oops:")


#==============================================================================
# End of extension
#==============================================================================
