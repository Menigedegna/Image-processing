# Imaris plugins for high resolution 3D nucleus microscopy images 

These plugins are designed for isolated nuclei undergone FISH or Immunostaining and conterstained with DNA staining dyes. The 3D microscopy images need to be high resolution for better performance (values in parameter files for each plugin can be adjusted to accomodate different resolution) and you should be able to visualise detect immunostaining foci or FISH spots visually on the images for the segmentation in the plugins to work. 
Please read the FISH_Image processing_protocol.pdf that I published in Methods and Protocols, Methods in Molecular Biology 2018, to get more information about how to acquire such (FISH) images.

## Setting up the  environment
Please follow the visual instructions in the Set_up.ppt to set up and use the plugins, and to have a better understanding on how the segmentation works, so that it is easier to adjust the input parameters for each plugin.

### Prerequisites

Imaris 8.4 and higher versions. 

Python 2.7 and higher versions and the following packages:
```
Tkinter
os
logging
tkMessageBox
numpy
pandas
tkFileDialog
datetime
shutil 

```

## Running the tests and getting Started

* You can download a demo file from the following link: 
```
https://www.dropbox.com/sh/7admke13knqzq3e/AADGlJVOT5_IZFcZLYTkfWffa?dl=0

```
Please run the plugins on the demo file provided to make sure that you have all the packages and software versions required to run the plugins properly. <br /> 
* Once you have you've confirmed that you get the desired result with the demo files, test few (10 images recommended) of your own images and adjust the parameters untill you get the desired results. <br /> 
* Once you have set and saved the input parameters, you can run an automated batch processing of all your images. Depending on the number and size of your images this may take a while. When you run any plugin a python console is displayed and keeps track of tasks completed and the duration. This can give you an idea of how much time is required for all your images.
* If you run into any problem or you have suggestions, please get in touch with me : mariamawitashenafi@yahoo.fr

## Built With

* Python 2.7
* Imaris 8.4 and tested on Imaris 9.2

## License

Please feel free to download, use and modify the plugins, during any publication, please attribute to Mariamawit S. Ashenafi. Thank you in advance and good luck.

## Acknowledgments

* UZH images processing plateform for providing Imaris software
* Dr Célia Baroux for collaborating with the nuclei isolation and FISH protocol
* Prof. Reinhard Furrer and Dr. Peter Majer for their recommendations
* Imaris developpers for the well organised and extensive documention on Imaris' build in functions.
* Users of Stack Overflow for their questions and usefull solutions.
