# Imaris plugins for high resolution 3D nucleus microscopy images 

These plugins are designed for isolated nuclei undergone FISH or Immunostaining and conterstained with DNA staining dyes. The 3D microscopy images need to be high resolution for better performance (values in parameter files for each plugin can be adjusted to accomodate different resolution) and you should be able to visualise detect immunostaining foci or FISH spots visually on the images for the segmentation in the plugins to work. 
Please read the FISH_Image processing_protocol.pdf that I published in Methods and Protocols, Methods in Molecular Biology, to get more information about how to acquire such (FISH) images.

## Setting up the  environment
Please follow the visual instructions in the Set_up.ppt to set up and use the plugins, and to have a better understanding on how the segmentation works, so that it is easier to adjust the input parameters for each plugin.

### Prerequisites

Imaris version 8.4 and above. 

Python 2.7 and above and the following packages:
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

Please run the plugins on the demo file provided to make sure that you have all the packages and software versions required to run the plugins properly. Once you have you confirm you get the desired result with the demo files, test few (10 images recommended) of your own images and adjust the parameters untill you get the desired results. Once you have set and saved the input parameters, run a batch processing of all your images, depending on the number and size of your images this may take a while. The python console keeps track of milestones achieved. This can give you an idea of how much time is required for all your images.

## Built With

* Python 2.7
* Imaris 8.4 and tested on Imaris 9.2

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* UZH images processing plateform for providing Imaris software
* Prof. Reinhard Furrer and Dr. Peter Majer for their recommendations
* Users of Stack Overflow for their questions and their answers
