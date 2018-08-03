# ==============================================================================
#
#    <CustomTools>
#      <Menu>
#       <Item name="XTSegment_nucleus" icon="Python" tooltip="XTSegment_nucleus">
#         <Command>PythonXT::XTSegment_nucleus(%i)</Command>
#       </Item>
#      </Menu>
#    </CustomTools>
# ==============================================================================

import ImarisLib
from Tkinter import *
import os
import logging
import tkMessageBox
import numpy as np
import pandas as pd
import tkFileDialog
import datetime
from shutil import rmtree

class Checkbar(Frame):
    """
    This is the Checkbar class, it sets up the pop up window
    """

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


class User_Set_Parameters():
    """
    This is the User_Set_Paraneters class, it creates a pop up window,
    and registers user's input.
    """

    def allstates(self):
        User_selection = list(self.lng.state())
        if sum(User_selection) > 0:
            self.root.destroy()
            self.SelectedOption = User_selection
        else:
            Message = "Please select one of the options."
            Label(self.root, text=Message).grid(row=1)

    def start_pop_window(self, Options=[], Message="No message"):
        self.root = Tk()
        Label(self.root, text=Message).grid(row=0)
        self.lng = Checkbar(self.root, Options)
        self.lng.grid(row=2)
        self.lng.config(relief=GROOVE, bd=2)
        Button(self.root, text='Quit', fg="red", command=quit).grid(row=4)
        Button(self.root, text='Submit', fg="darkgreen", command=self.allstates).grid(row=5)
        self.root.mainloop()

    def __init__(self, selection_options=[], message="None"):
        self.start_pop_window(selection_options, message)


class Imaris_Plugin():
    """
    This is the ImarisApp class, it connects to an imaris instance,
    registers static data : plugin name and parameters required by parameters,
    and keeps track of processed files
    """
    vImaris = None
    file_name_list = None

    def error_in_csv_file(self):
        tkMessageBox.showinfo(title="Alert",
        message="Please make sure that the " + self.plugin_name + "_Parameters.csv file contains the following parameters in this order:"
        + "\n- Surface grain value for detailed nucleus segmentation : float"
        + "\n- Surface grain value for rough nucleus segmentation : float"
        + "\n- Surface grain value for nucleolus segmentation : float"
        + "\n- Local surface diameter for nucleolus segmentation : float"
        + "\n- Local surface diameter for nucleus segmentation : float"
        + "\n- Surface grain value for chromocenter segmentation : float"
        + "\n- Local surface diameter for chromocenter segmentation : float"
        + "\n- Remove nucleus surface close to border : boolean"
        + "\n- Select nucleus surface with highest volume : boolean"
        + "\n- Minimum volume threshold for chromocenter segmentation: float"
        + "\n- Minimum intensity threshold for chromocenter segmentation: 0 < float < 1"
        + "\n- Minimum volume threshold for nucleus segmentation: float"
        + "\n- Save image : boolean"
        + "\n- Sigma value for gaussian filter channel smoothing"
        + "\n- Sigma value for channel background substraction")
        quit()

    def get_plugin_parameters(self):
        currentDirectory = os.getcwd()
        AllFilesInDirectory = os.listdir(currentDirectory)
        ParametersList = None
        parameter_file = self.plugin_name+"_Parameters.csv"
        if parameter_file in AllFilesInDirectory:
            ParameterData = pd.read_csv(parameter_file, sep=";", header='infer', decimal='.')
            if "Value" in ParameterData.columns:
                ParametersList = list(ParameterData["Value"])
                return ParametersList
            else:
                tkMessageBox.showinfo(title="Error", message="Please make sure the '" + parameter_file + "' file contains a column 'Value' containing the values necessary for this plugin.")
                quit()
        else:
            tkMessageBox.showinfo(title="Error", message="Please make sure there is a '" + parameter_file + "' in the folder containing this plugin.")
            quit()

    def get_saved_parameters(self):
        parameter_list = self.get_plugin_parameters()
        parameter_list = [float(x) for x in parameter_list]
        if len(parameter_list) >= 15:
            self.nucleus_detail_surface_grain = parameter_list[0]
            self.nucleus_rough_surface_grain = parameter_list[1]
            self.nucleolus_surface_grain = parameter_list[2]
            self.nucleolus_surface_surface_diameter = parameter_list[3]
            self.nucleus_surface_surface_diameter = parameter_list[4]
            self.chromocenter_surface_grain = parameter_list[5]
            self.chromocenter_surface_diameter = parameter_list[6]
            self.remove_border_toggle = parameter_list[7]
            self.select_biggest_nucleus_toggle = parameter_list[8]
            self.chromocenter_minimum_volume = parameter_list[9]
            self.chromocenter_minimum_intensity = parameter_list[10]
            self.nucleus_minimum_volume = parameter_list[11]
            self.save_image_toggle = parameter_list[12]
            self.gaussian_filter_sigma = parameter_list[13]
            self.background_substraction_sigma = parameter_list[14]
        else:
            self.error_in_csv_file()

    def get_user_set_parameters(self):
        pop_up_instance = User_Set_Parameters(
        ["Batch of images", "Just one image"],
        "Do you wish to run the script on a batch of images or just on opened image?"
        )
        # get first boolean: is "Batch of images" selected?
        self.BatchProcessingToggle = pop_up_instance.SelectedOption[0]

        pop_up_instance = User_Set_Parameters(
        ["Yes", "No"],
        "Would you like to run the segmentation module?"
        )
        # get first boolean : is "yes" selected?
        self.SegmentationToggle = pop_up_instance.SelectedOption[0]

        pop_up_instance = User_Set_Parameters(
        ["The nucleus", "The nucleolus", "The chromocenters"],
        "Which surface(s) are you analysing?"
        )
        # get list of boolean
        self.SurfaceOption = pop_up_instance.SelectedOption

    def get_user_set_directory(self):
        root1 = Tk()
        Image_folder = tkFileDialog.askdirectory(parent=root1, initialdir="/", title='Please select the directory containing the images to be processed. \n The folder containing the resulting files will be saved in this directory.')
        root1.destroy()
        self.image_folder = Image_folder

    def get_file_directory(self):
        vFileName = self.vImaris.GetCurrentFileName()
        vFilePath = os.path.dirname(vFileName)
        self.image_folder = vFilePath
        return os.path.split(vFileName)[1]

    def create_directory(self):
        Result_pathway = os.path.join(self.image_folder, self.plugin_name + "_Result")
        counter = 0
        while os.path.exists(Result_pathway) and counter < 3:
            # tkMessageBox.showinfo(title="Alert",
            # message="Please save the folder '" + self.plugin_name + "_Result' under another name first!")
            # quit()
            print("Please save the folder '" + self.plugin_name + "_Result' under another name!")
            print("If you press ENTER three times, without renaming your folder, IT WILL BE DELETED!")
            raw_input("Press ENTER to continue.")
            counter += 1
        if os.path.exists(Result_pathway):
            rmtree(Result_pathway)
        os.makedirs(Result_pathway)
        self.result_path = Result_pathway


    def __init__(self, aImarisId, plugin_name):
        try:
            vImarisLib = ImarisLib.ImarisLib()
            self.vImaris = vImarisLib.GetApplication(aImarisId)
            self.plugin_name = plugin_name
            if self.vImaris is not None:
                # Get parameters saved in .csv file with the same name as this plugin
                self.get_saved_parameters()
                # Set up pop up windows to ask user for few parameters
                self.get_user_set_parameters()
                if self.BatchProcessingToggle:
                    self.get_user_set_directory()
                    self.file_name_list = [i for i in os.listdir(self.image_folder) if i.endswith('.ims') or i.endswith('.ics')]
                else:
                    vFileName = self.get_file_directory()
                    self.file_name_list = [vFileName]
                self.create_directory()
            else:
                tkMessageBox.showinfo(title="Alert", message="Can't connect with Imaris!")
                quit()
        except:
            logging.exception("Oops: Error inside class Imaris_Plugin")


class Images():

    _plugin_parameters = None
    _number_processed_file = 0

    def get_plugin_instance(self, plugin_instance):
        if self._plugin_parameters is None:
            self._plugin_parameters = plugin_instance

    def get_channel_number(self):
        self.Total_channel = self.image.GetSizeC()

    def get_selected_channel(self):
        self.get_channel_number()
        pop_up_instance = User_Set_Parameters(
        range(self.Total_channel),
        "Select DAPI channel.")
        # get first occurance of True in list
        dapi_channel = [i for i, x in enumerate(pop_up_instance.SelectedOption) if x == 1]
        self.dapi_channel = dapi_channel[0]

    def logtime(self, task_name="", give_duration=True):
        # This is the logtime, logs start and end of tasks
        curtime = datetime.datetime.now()
        if self.gLasttime is not None and give_duration:
            diff = (curtime-self.gLasttime).total_seconds()
        else:
            diff = curtime.ctime()
            self.gLasttime = curtime
        print(str(diff) + ' sec. : ' + task_name)

    def open_file(self):
        full_path_file_name = os.path.join(self._plugin_parameters.image_folder, self.file_name)
        self._plugin_parameters.vImaris.FileOpen(full_path_file_name, "")


    def __init__(self, plugin_instance, file_name):
        self.gLasttime = None
        self.file_name = file_name
        try:
            self.get_plugin_instance(plugin_instance)
            self.open_file()
            self.image = self._plugin_parameters.vImaris.GetDataSet()
            if self.image is not None:
                task_name = "Start processing file - " + self.file_name
                self.logtime(task_name=task_name)
                self.volume = self._plugin_parameters.vImaris.GetSurpassSelection()
                self.factory = self._plugin_parameters.vImaris.GetFactory()
                self.scene = self._plugin_parameters.vImaris.GetSurpassScene()
                self.container = self.factory.CreateDataContainer()
                self.container.SetName('Segmented objects')
                if self._number_processed_file == 0:
                    self.get_selected_channel()
            else:
                print("No image detected in file" + self.file_name)
        except:
            logging.exception("Oops: Error inside class Images")

    def remove_container(self):
        self.scene.RemoveChild(self.container)

    def exit_file(self):
        if self._plugin_parameters.save_image_toggle:
            full_path = os.path.join(self._plugin_parameters.result_path, self.file_name+".ims")
            self._plugin_parameters.vImaris.FileSave(full_path, "")
            task_name = "Save file"
            self.logtime(task_name=task_name)
        self.remove_container()
        task_name = "End processing file - " + self.file_name
        self.logtime(task_name=task_name, give_duration=False)
        self._number_processed_file += 1

# ==============================================================================
#  FUNCTIONS REQUIRED BY SURFACE CLASS
# ==============================================================================
# Get numeric data for an item in scene
def GetStat(vSceneItem, FilterString):
    vAllStatistics = vSceneItem.GetStatistics()
    vNames = vAllStatistics.mNames
    vValues = vAllStatistics.mValues
    OutputParam = [float(vValues[a]) for a, x in enumerate(vNames) if x == FilterString]
    return OutputParam

# Get max intensity from list of intensity values #
def getMax(ListIntensity):
    maxList = [max(i) for j in ListIntensity for i in j]
    maxValue = max(maxList)
    return maxValue

# Select voxels with low DAPI intensity #
def getLowIntensity(detailed_nucleus_mask_values, smoothed_channel_mask, rough_nucleus_mask_values):
    maxValue = getMax(smoothed_channel_data)
    Result = []
    for z in range(len(smoothed_channel_data)):
        mask_y = []
        for y in range(len(smoothed_channel_data[0])):
            mask_surface = detailed_nucleus_mask_values[z][y]
            mask_surface2 = rough_nucleus_mask_values[z][y]
            mask = smoothed_channel_mask[z][y]
            # DAPI MASK: NORMALISE INTENSITY WITH MAXIMUM VALUE AND ROUND UP TO 1
            mask = [round(item / maxValue, 1) for item in mask]
            # = SELECT VOXELS WITH BELOW AVERAGE INTENSITY
            mask = [0 if item >= 0.5 else item for item in mask]
            # = SELECT VOXELS OUTSIDE OF DETAILED SURFACE
            mask_surface = [maxValue if item == 0 else 0 for item in mask_surface]
            # REMOVE VOXELS OUTSIDE OF ROUGH SURFACE
            mask = [i * j * k for i, j, k in zip(mask, mask_surface, mask_surface2)]
            mask_y.append(mask)
        Result.append(mask_y)
    return Result


class Surface(Images):
    """
    This is Surface class. It inherits Images class and segments the nucleus,
    the nucleolus and the chromocenters
    """

    def create_surface_in_scene(self):
        # ADD ITEM TO SCENE
        self.surface.SetName(self.object_type)
        self.container.AddChild(self.surface, -1)
        self.scene.AddChild(self.container, -1)
        # SAVE SNAPSHOT AND DESELECT ITEM FROM SCENE
        full_path = os.path.join(self._plugin_parameters.result_path,
        "Snapshot" + self.object_type +"_"+self.file_name+".tif")
        self._plugin_parameters.vImaris.SaveSnapShot(full_path)
        #deselect item from scene
        self.surface.SetVisible(0)

    def select_surface(self):
        vTimeIndex = 0
        vol = GetStat(self.surface, "Volume")
        if self.object_type == "nucleolus":
            area = GetStat(self.surface, "Area")
            volumeToAreaRatio = [x/y for x, y in zip(vol, area)]
            SelectedID = max(xrange(len(volumeToAreaRatio)), key=volumeToAreaRatio.__getitem__)
        if self.object_type == "nucleus":
            SelectedID = max(xrange(len(vol)), key=vol.__getitem__)
        vertices = self.surface.GetVertices(SelectedID)
        vNormals = self.surface.GetNormals(SelectedID)
        faces = self.surface.GetTriangles(SelectedID)
        self.surface = self.factory.CreateSurfaces()
        self.surface.AddSurface(vertices, faces, vNormals, vTimeIndex)

    def set_string_filter(self):
        SfS = '"Volume" above ' + str(self._plugin_parameters.nucleus_minimum_volume) + ' um^3'
        if self._plugin_parameters.remove_border_toggle:
            SfS = SfS + ' "Distance to Image Border XYZ" above 0.0516 um'
        return SfS

    def segment_surface(self, surface_grain, channel_index, string_filter, local_diameter, vATM):
        vATA = 1
        if vATM > 0:
            vATA = 0
        vROI = None
        surface = self._plugin_parameters.vImaris.GetImageProcessing().DetectSurfaces(
                    self.image,
                    vROI,
                    channel_index,
                    surface_grain,
                    local_diameter,
                    vATA,
                    vATM,
                    string_filter)
        return surface

    def Get_Mask_data(surface):
        vImageSizeX = self.image.GetSizeX()
        vImageSizeY = self.image.GetSizeY()
        vImageSizeZ = self.image.GetSizeZ()
        vExtentMinX = self.image.GetExtendMinX()
        vExtentMinY = self.image.GetExtendMinY()
        vExtentMinZ = self.image.GetExtendMinZ()
        vExtentMaxX = self.image.GetExtendMaxX()
        vExtentMaxY = self.image.GetExtendMaxY()
        vExtentMaxZ = self.image.GetExtendMaxZ()
        mask_min = [vExtentMinX, vExtentMinY, vExtentMinZ]
        mask_max = [vExtentMaxX, vExtentMaxY, vExtentMaxZ]
        mask_size = [vImageSizeX, vImageSizeY, vImageSizeZ]
        mask_time =	0
        mask =	surface.GetMask(mask_min[0], mask_min[1], mask_min[2], mask_max[0], mask_max[1], mask_max[2],mask_size[0],mask_size[1], mask_size[2], mask_time)
        mask_values =	mask.GetDataVolumeFloats(0,0)
        return mask_values

    def AddChannel(self, list_voxel_intensity, channel_name):
        time_index = 0
        channel_index = self.Total_channel
        DAPI_channel_color = self.image.GetChannelColorRGBA(self.dapi_channel)
        self.Total_channel += 1
        self.image.SetSizeC(self.Total_channel)
        self.image.SetDataVolumeFloats(list_voxel_intensity, channel_index, time_index)
        self.image.SetChannelName(channel_index, channel_name)
        self.image.SetChannelColorRGBA (channel_index, DAPI_channel_color)

    def new_segment_decorator(segment_func):
        def edit_scene(self):
            segment_func(self)
            if self.surface is not None:
                if self.object_type in ["nucleus", "nucleolus"]:
                    self.select_surface()
                self.create_surface_in_scene()
            else:
                print(self.object_type.upper() +" IS NOT DETECTED!")
        return edit_scene

    @new_segment_decorator
    def detect_nucleolus(self):
        surface = None
        vATM = 0
        string_filter=self.set_string_filter()
        # Do a detailed nucleus segmentation: set a low sigma to gaussian filter
        detailed_nucleus_segmentation = self.segment_surface(
        self._plugin_parameters.nucleus_detail_surface_grain,
        self.smooth_dapi_channel,
        string_filter,
        self.nucleus_surface_surface_diameter,
        vATM)
        if detailed_nucleus_segmentation.GetNumberOfSurfaces() > 0:
            # Do rough nucleus segmentation: set a high sigma to gaussian filter
            rough_nucleus_segmentation = self.segment_surface(
            self._plugin_parameters.nucleus_rough_surface_grain,
            self.smooth_dapi_channel,
            string_filter,
            self._plugin_parameters.nucleus_surface_surface_diameter,
            vATM)
            # Get DAPI channel mask
            smoothed_channel_data = self.image.GetDataVolumeFloats(self.smooth_dapi_channel,0)
            # Get mask for detailed / rough nucleus segmentation
            Detailed_nucleus_mask = Get_Mask_data(detailed_nucleus_segmentation)
            Rough_nucleus_mask = Get_Mask_data(rough_nucleus_segmentation)
            # From the above masks, infere nucleolus mask
            Nucleolus_mask = getLowIntensity(
            Detailed_nucleus_mask,
            smoothed_channel_data,
            Rough_nucleus_mask)
            # Create new channel with nucleolus mask
            self.AddChannel(Nucleolus_mask, "Low DAPI intensity")
            # Segment new channel
            string_filter = ""
            surface = self.segment_surface(
            self._plugin_parameters.nucleolus_surface_grain,
            self.Total_channel-1,
            string_filter,
            self._plugin_parameters.nucleolus_surface_surface_diameter,
            vATM)
        self.surface = surface

    @new_segment_decorator
    def detect_chromocenters(self):
        numberSceneInstance = self.scene.GetNumberOfChildren()
        VolumeNotFound = True
        i = 0
        while i <= numberSceneInstance and VolumeNotFound:
            selection = self.scene.GetChild(i)
            aVolume = self.factory.ToVolume(selection)
            if aVolume is not None:
                VolumeNotFound = False
            i += 1
        IntensityMax = GetStat(aVolume, "Data Intensity Max")[self.smooth_dapi_channel]
        string_filter = '"Volume" above '+str(self._plugin_parameters.chromocenter_minimum_volume)+' um^3'
        vATM = IntensityMax * self._plugin_parameters.chromocenter_minimum_intensity
        surface_grain = self._plugin_parameters.chromocenter_surface_grain
        local_diameter = self._plugin_parameters.chromocenter_surface_diameter
        surface = self.segment_surface(
        surface_grain,
        self.smooth_dapi_channel,
        string_filter,
        local_diameter,
        vATM)
        self.surface = surface if surface.GetNumberOfSurfaces() > 0 else None

    @new_segment_decorator
    def detect_nucleus(self):
        vATM = 0
        string_filter = self.set_string_filter()
        surface = self.segment_surface(
        self._plugin_parameters.nucleus_rough_surface_grain,
        self.dapi_channel,
        string_filter,
        self._plugin_parameters.nucleus_surface_surface_diameter,
        vATM)
        self.surface = surface if surface.GetNumberOfSurfaces() > 0 else None

    def apply_gaussian_filter(self):
        vImage_data = self.image.GetDataVolumeFloats(
        self.dapi_channel,0)
        self.AddChannel(vImage_data, "Smooth DAPI channel")
        self.smooth_dapi_channel = self.Total_channel-1
        self._plugin_parameters.vImaris.GetImageProcessing().GaussFilterChannel(
        self.image,
        self.smooth_dapi_channel,
        self._plugin_parameters.gaussian_filter_sigma)
        self._plugin_parameters.vImaris.GetImageProcessing().SubtractBackgroundChannel(
        self.image,
        self.smooth_dapi_channel,
        self._plugin_parameters.background_substraction_sigma)

    # TODO: GET FEATURES OF SURFACES
    # TODO: CHECK NUCLEOLUS AND CC DETECTION
    # TODO: ERROR: END FILE PROCESSING DISPLAYED BEFORE END OF CLASS SURFACE

    def Segment_nucleus(self):
        try:
            # Segment the nucleus
            self.object_type = "nucleus"
            self.detect_nucleus()
            self.nucleus = self.surface
            self.logtime(task_name="Detect nucleus surface")
            if self._plugin_parameters.SurfaceOption[1] or self._plugin_parameters.SurfaceOption[2]:
                # Apply gaussian filter to DAPI channel and create new channel
                self.apply_gaussian_filter()
                self.logtime(task_name="Smooth DAPI channel")
            if self._plugin_parameters.SurfaceOption[1]:
                # Segment nucleous
                self.object_type = "nucleolus"
                self.detect_nucleolus()
                self.logtime(task_name="Detect nucleus surface")
            if self._plugin_parameters.SurfaceOption[2]:
                # Segment chromocenters
                self.object_type = "chromocenters"
                self.detect_chromocenters()
                self.logtime(task_name="Detect chromocenter surface")
        except:
            logging.exception("Oops: Error inside class Surface")

def XTSegment_nucleus(aImarisId):
    # Connect to Imaris
    plugin_name = "XTSegment_nucleus"
    logging.basicConfig(level=logging.DEBUG, filename="log["+plugin_name+"].log")
    try:
        plugin_instance = Imaris_Plugin(aImarisId, plugin_name)
        for file_name in plugin_instance.file_name_list:
            # image_instance = Images(plugin_instance, file_name)
            # print(image_instance._number_processed_file)
            image_instance=Surface(plugin_instance=plugin_instance, file_name=plugin_instance.file_name_list[0])
            image_instance.Segment_nucleus()
            print("plugin working")
            raw_input("Press Enter to terminate.")
    except:
        logging.exception("Oops: Error in main XT function")
