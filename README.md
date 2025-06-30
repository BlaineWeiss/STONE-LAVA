# STONE / LAVA - Scientific Analysis Software

Copyright © 2025 Blaine Everett Weiss, University of Kentucky Research Foundation
All rights reserved.

This software is distributed for non-commercial academic research use only.
Use in any commercial setting, including for-profit institutions or fee-for-service labs,
is strictly prohibited without a separate commercial license.

A patent has been filed on the underlying method. Commercial use may infringe on that IP.
Patent application number: Pending

License terms: See LICENSE_ACADEMIC.txt and LICENSE_COMMERCIAL.txt for details.

For commercial licensing inquiries, contact: blaine.weiss@uky.edu

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

If you use this software in published work, please cite:

Disrupted calcium dynamics in reactive astrocytes occur with endfeet-arteriole decoupling in an amyloid mouse model of Alzheimer’s disease
Blaine E. Weiss, et al. bioRxiv 2025.01.24.634584; doi: https://doi.org/10.1101/2025.01.24.634584

STONE-LAVA is a combined analysis platform for the quantitation of functional imaging data acquired by two-photon microscopy. This document will provide a general overview of their key functions, and their utilization for comparision of cell activity to vascular dynamics.  

#### STONE: (Spatial & Temporal Observation of Network Events) aids in the identification and subcompartmentalization of active cells. Fluorescence time series data extracted from these segmentations are analyzed for signaling events, and transient parameters such as amplitude, rise/decay times, and area under the curve (AUC) are determined. Once all events are indexed, synchronicity of cells within the network is determined by those events and are mapped.

#### LAVA: (Localized Analysis of Vascular Astroyctes) is a standalone application with optional pairing with STONE acquired data. LAVA's initial function is to quantify vascular tone changes during imaged experiment trials. It does so by modeling vascular architecture and then measuring cross sectional diameters of vessels of interest. Next, the results can be used to quantify fluorescence changes that occur in perivascular spaces immediately adjacent to the vessel of interest. This allows for directional and correlational assessment of vascular tone changes and vascular/perivascular cell activity. When launched from STONE, cell segmentation data may be used to make additional comparisions with vascular and perivascular activity from the same trial, and field of view.

## Citation/Contact Information
We sincerely hope you will find value in these applications in which case we ask that you please include the following reference in your publications:

Weiss et al.......

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">Spatial & Temporal Observation of Network Events (STONE) & Localized Analysis of Vascular Astrocytes (LAVA) software </span> is licensed under <a href="https://creativecommons.org/licenses/by-nc-sa/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY-NC-SA 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/nc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/sa.svg?ref=chooser-v1" alt=""></a></p>

We are continually working to improve the performance and versatility of the techniques used in this package. Therefore we encourage users to submit recommended improvements through commits. Pull requests will be considered based on compatiability and benchmark improvements for a particular technique.

For specific inquiries or recommended improvements to these applications, also feel free to contact us at
Blaine.Weiss@uky.edu

## Table of Contents
1. [Hardware Requirements & Optimization Information](#Hardware-Requirements--Optimization)
2. [STONE Layout](#STONE-Layout)
3. [Loading/Importing imaging and .mat data](#LoadingImporting-imaging-and-mat-data)
4. [Preprocessing/Analysis selection](#PreprocessingAnalysis-selection)
    - Motion Correction
    - Downscaling
    - Spatial Filtering
5. [Detecting regions of interest](#Detecting-regions-of-interest)
6. [Troubleshooting and help](#Troubleshooting-and-help)
7. [Event Detection & Signal Visualization](#Event-Detection--Signal-Visualization)
8. [Event Editor & Parameter Collection](#Event-Editor--Parameter-Collection)
9. [Network Synchronicity & Frequency Analysis](#Network-Synchronicity--Frequency-Analysis)
10. [Saving and Data Structuring](#Saving-and-Data-Structuring)
11. [MATLAB Workspace Integration](#MATLAB-Workspace-Integration)
12. [LAVA](#LAVA)
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>


## Hardware Requirements & Optimization
Performance optimizations are currently being conducted to improve computational efficiency, and lower the required hardware specifications.


Analysis machines will require a m
Minimum recommended specs: 

Processor: 8 core Intel Xeon(R) W-2145 CPU or greater
GPU: NVIDIA Quadro P1000
RAM: 32GB for files of size 2.5 GB


## STONE Layout
<a name= "AppLayout"></a>
Custom link
<br/>
<br/>
<br/>

![image](https://github.com/user-attachments/assets/f885aa8e-23de-404f-996c-3f3b62ee4ebd)

1 File name display\
2 Start motion correction function\
3 Downscaling bin factor\
4 Initiate/refresh image stack with specified downscale factor\
5 Application status label\
6 Application busy indicator\
7 Flattened activity scaled image\
8 Image Calibration [See Metadata](#Metadata)\
9 Application Log\
10 Segmentation Label\
11 Make Filter Command\
12\
1\

1 


First, the typical workflow for STONE will be described, with alternative workflows and options described thereafter.

## Loading/Importing imaging and .mat data
### Load File
Clicking the "Load File" button will call a dialog window to select a file to analyze. .Tif, .Tiff, and .mat are currently the only accepted formats. Current import format is a 3D interleaved channel matrix following M x N x S order. 

**4D array formatting is in progress

After array import, the user is prompted to select a channel which represents the cell fluorescence activity of interest. Once selected, all channels will be analyzed and corrected for bidirectional scanning offsets that may exist.

Existing .mat data files exported from prior application instances can also be loaded through this function to import previous segmentation data, and/or metadata. See [Saving and Data Structuring](#Saving-and-Data-Structuring) for more details. After .mat files are loaded, a prompt for the assignment of data to a segmentation label will appear, and masks will be stored in the application instance. You may repeat the Load function until all desired segmentations are loaded. The user may then proceed to additional analysis. Please note that prior analysis of somata segmentations is generally recommended before fine process segmentation. See [Detecting regions of interest](#Detecting-regions-of-interest) for details. 



### Metadata
 Tag based metadata collected during acquisition of .tif experiment files can be loaded into the application workspace using the "Load File" command. Relevant information such as image scale, acquisition frame rate, and frame dimensions are required for size and signal parameter measurement. Currently metadata obtained from ScanImage .tif files will automatically populate the necessary variables for calibration. In the absence of tag metadata, the application will prompt for the entry of calibration information.

 For specific tag structure inclusion, you can edit the Load File function variables to match your tag fields, or [contact us](#CitationContact-Information) and include details of your imaging software platform, and an example data file containing metadata.
 
 
 
 <a% name="Metadata"></a>

## Preprocessing/Analysis selection
### Motion Correction
For awake animal imaging, it is sometimes necessary to perform motion correction to stabilize the field of view for accurate ROI identification. The application utilizes the commonly used NormCorre algorthim for image registration and stabilization.


### Downscaling
For activity dependent identication of active cell features, the application will extract time series fluorescence data from all pixel points for a 2D FOV. For data size and noise reduction purposes, it is beneficial to downscale the data array in the M & N dimensions. The blocksize value determines the single dimension averaging bin length in current pixel values. Downscaling will always occur equally in M & N dimensions, and assumes square pixel sizes from acquisition. Depending on native resolution, pixel binning may result in a small remainder of ungrouped pixels on the image edges. Those edges are trimmed off the array and image calibration is adjusted for the cropped array, and represented at the bottom of the image display ([8](#Stone-Layout))

%Diagram of 3x3 block 4x4 and plotted exponential data reduction.

Increasing the block size as shown in panel A, results in an exponential decrease in signal array size. (Panel B) For a 512x512 FOV a factor of 3 will bin pixels into 3x3 groups, resulting in a 170x170 FOV with 1 pixel removed from all edges. Time series values for each pixel bin are the average of the 9 pixel values that make up the bin.


###


## Detecting regions of interest
### Image Filtering
Depending on magnification and chosen cell label, discerning cell boundaries for ROI detection can be challenging, particularly for activity dependent identification using genetically encoded Ca2+ indicators. Astrocyte connectivity via gap junctions cause signals to be synchronuous and propagative in nature, and can significantly confound measures of event kinetics. STONE's ROI detection performance is significantly augmented by the application of digital filters to the application array. Below is a detailed description of the technique and instructions for use.

Filter creation tool is a standalone application that creates custom contrast enhancing and size attenuating/accentuating filters. When used with the STONE application, these filters non-destructively separate array features desired for signal analysis. The application can be called directly from the STONE GUI "Make Filter" button [11](#Stone-Layout) The filter tool will load with the designated fluorescence channel as its primary input. Once the representative activity image is loaded, the user may begin using the filter application. 


%%% Add visual of app
%%% Label image
%%% Add instructions for use.

The "Load Filter" button on the bottom right of the application may be used to load preset filters downloaded with the package. The "Soma Filter" is ideal for separation of closely connected cells. Once loaded, the filter can be tuned by checking the "Smooth Filter" checkbox and adjusting the slider below. Crucially in trials with unlevel Z planes, the "Sharp Center" command should be used to effectively denoise parenchymal space. Core filter parameters such as band pass can be adjusted above, but will require reprocessing of the filter function. This is completed by pressing the "Form Filter" button. More details on filter creation and modification are explained below.


How to make a custom filter: (See Filter Creation README)


### GUI Interactions with Filter Creation Tool:

To push the filter to the main application, press the "Send Filter" button. The main app status label will read "Filter Received" to confirm that the data was transferred successfully. To clear the filter you may press "Clear Filter" on the main application. The filter used by the main applicaton can be replaced by an updated one from the filter creation application by pressing the "Send Filter" button again. 

Between parent application instances, the filter creation app may remain open for convenience if the user is using the same filter over multiple data files. The filter tool will by default send filters and app data to the parent "STONE" application.

This is further explained below in [Voxel Clustering](#Voxel-Clustering)


### Voxel Clustering
Once all preprocessing and filtering steps have been completed, voxel clustering may be performed for segmentation. The "Detect ROIs" button will initiate the clustering algorithm. Checking the "Intensity Filter" checkbox above will include an intensity thresholding step prior to clustering. This is useful if filtering was not effective in background subtraction, the file has bad signal to noise ratio, or if activity is small from cells of interest. The threshold is set to retain voxels 1 standard deviation above the mean.
%%%confirm threshold

Once the "Detect ROIs" button is pressed, raw signal and signal parameters from each voxel are calculated and input into a principal component analysis. The top 4 components are then clustered by kmeans to output 3 proposed segmentations. The result is clear separation between active and non active pixels in the FOV. 2/3 cluster groups represent regions of cell activity. Moving along the largest principal component typically results in voxels clustering further away from cell somas. For best accuracy, it is recommended to select the cluster group containing cell somata first. This is because subsequent fine process segmentation will mask away all prior segmentations, including somata, and result in more sensitive detection of fine process signaling.


Here include the clustering results...





![ClusterSelection](https://github.com/user-attachments/assets/302c8fc1-aef5-491f-a499-134338c44e8e)

Scaled Activity Image        |    Somata Segmentation Result
:----:|:--:
![ScaledActivityImage208spon2_5](https://github.com/user-attachments/assets/1bfed750-4ef9-4ad7-9054-decd5a7ae070)    |    ![ROIsSOMA](https://github.com/user-attachments/assets/cfc32bf9-ef7f-4bb4-b5bf-15b9b8ac4fce)


Watershed splitting

![ROIsSOMASplit_h_adjustedWatershed](https://github.com/user-attachments/assets/cbdc4f2b-6760-494f-a6bc-83ab281dc305)


A file example to select for cell somata.



### Split ROI function

After clustering is completed, the resulting ROIs can be further refined by separating contiguous cells that may have been clustered as one. Pressing the "Split ROIs" button and selecting a ROI enveloping multiple cells will first attempt separation using the watershed transform. If the ROI fails to split, manual splitting is enabled. To manually split a ROI, select click on a point outside or along the perimeter of the ROI where the split should occur. Then select draw points that define the new boundary. The tool will assist in selecting endpoints on the perimeters of the ROI. Once the boundary is completed, right click to finalize the split path. If you make a mistake, you may press the "Undo" button next to the "Split ROIs" button to revert the ROI to its previous shape.

place image or video gif


### Perivascular Activity Segmentation
If imaging data includes a channel containing vasculature, you may restrict the segmentation algorithm to cluster only the pixels on or around the perivascular space. This is useful for the identication of astrocyte endfeet, pericytes, and other various cell types localized to the neurovascular unit. Changing the Analysis Type to "Perivascular" will launch a vascular thresholding tool that will assist in masking image stacks to perivascular space. Upon tool opening, a flattened vascular channel image will pass to the display.

Using the "Thresholding" slider, an intensity threshold can be applied to mask the vascular space. Depending on imaging location and conditions, an uneven imaging field may cause difficulty using global thresholds. In these cases it is recommended to enable the "Local Contrast Enhancement" for improved vessel detection.

%%%display angled FOV and local contrast enhancement








<table width="100%">
<thead>
<tr>
    <th width="33%">Channel Composite Image
    <th width="33%">Global Threshold Mask
    <th width="33%">Filtered & Thresholded Mask
<tr>
<tr>
    <td width="33%"><img src=https://github.com/user-attachments/assets/44dc1810-f0b5-40d4-9b9d-1b859706d454/></td>
    <td width="33%"><img src=https://github.com/user-attachments/assets/e9491283-178f-444c-8910-ecb284142a39/></td>
    <td width="33%"><img src=https://github.com/user-attachments/assets/a086120d-42f0-445d-a6b9-a3fa9a187913/></td>
</tr>
<thead>
<tbody>
    <td width="33%">Vascular channel is viewed by its flattened mean image
    <td width="33%">Global thresholding is applied by moving the "Thresholding" slider to include the desired regions. "Local Contrast Enhancement will restrict the degree of change by the slider.
    <td width="33%">If cross spectra contamination, imaging artifacts, or bad signal to noise prevents effective vascular masking, you can utilize the Image Filtering tool described (previously) to eliminate or attenuate undesired features. 
<tbody>
<table>





## Guidance for Complexity in Imaging FOVs
A common challenge in analysis involving segmentation is decisiveness on the handling of heterogeneous cell morphology and FOVs. In the case of astrocyte compartment segmentation, it may be difficult to decipher cell compartments for cells in close proximity to vascular spaces. While many astrocytes contain distal processes that terminate to endfeet, others may entirely ensheath vascular spaces, therby contaminating endfoot classification with cell soma.


Here is an example of a FOV enriched with vessel localized astrocytes. Somata and fine process segmentation can still be performed. Additional parsing of perivascular Ca2+ will necessarily have overlapping segmentation with somata ROIs.

Channel Composite Image             |  Segmentated Astrocyte Subcompartments
:-------------------------:|:-------------------------:
![LAVAdemo composite](https://github.com/user-attachments/assets/9349d661-a5d3-4001-87f5-df681a1d4da7)  | ![144file14c1somagliopilsegmentationsWhiteRed](https://github.com/user-attachments/assets/382bd773-6e2c-49ad-9ca4-dd592a40e469)


Depending on the research question, it may be necessary to compare subcompartment analysis parameters to other cell types/morphological structures. See LAVA ####LINK#### for comparisions to vascular dynamics and blood flow.

## Event Detection & Signal Visualization

Once regions of interest have been identified and indexed, time series vectors of pixels contained for each ROI are averaged together. The resulting time series vector represents the fluorescence intensity of the ROI over time.

A time index of ROI fluorescence activity can be used to quantify the signaling properties of the cells of interest, and to correlate interconnected brain processes, and signaling that indicate functional connectivity.
To acquire this index, another convolutional filtering algorithm is used on the ROI signal vectors to determine time regions of active cell activity.

Due to a potentially high number of active cells, the user has the option to restrict the number of ROIs used to save compute time and/or data size. The value of the numerical type box labeled "Maximum ROIs" can be edited directly or by using the up and down arrows attached to the box. Once the number of ROIs is determined, pressing "Trace ROI Extraction" will begin the event detection function.

### Trace ROI Extraction:

- Convolutional filtering technique and thresholding:
  
    For extraction of significant fluorescence events of variable amplitude and/or duration, a technique using wavelet transformation and convolutional filtering is applied on ROI signal extracts within an expected boundary of event timing. First, candidate events are determined by parallel convolution of the signal array with model sinusodal waveforms of various frequencies. The results are normalized and structured into matrix form. Next, the power spectrum computed by wavelet transformation is evaluated over these boundaries, smoothed to remove impulse artifacts and noise, and is then used to weigh candidate events in the time series arrays.

  $$ \sum_{i=1}^n wt(i,:)_{smoothed} * y(i) $$

Variable | Definition
:--:|:----:
  n | number of convolutions
  wt | Wavelet Transform Matrix
  y | Convolution Outputs
  k | Event frequency iterator

  
Negative shifts over 0.5 STD in the convolution are attenuated prior to weight combinations. After combination peak detection is performed with an initial screening threshold of 0.5 STD, the index start and end times are determined.









## Event Editor & Parameter Collection
Display

Consider repeating app layout image here...



## Network Synchronicity & Frequency Analysis

https://github.com/user-attachments/assets/0ebe409d-1bca-40c7-8fae-e0b08414e1a2




## Saving and Data Structuring
Clicking the Save button at the bottom of the application will open a command menu to save results, motion corrected image sequences, and reduced data .tifs for video production.

![image](https://github.com/user-attachments/assets/4503cca7-7b67-49cc-91b7-6eb5c5b9ae74)

### Output Data
The Output Data command will begin saving the results of all processed segmentations along with their event indices, signal parameters, and correlations. By default, the data will be saved in a nested directory named 

"ROI Segmentation Analysis\ _Segmentation Type_ \ _Filename_ .mat" 

This format allows for easy batch consolidation of trial results with segmentation analysis of the same type. Segmentation Type will be either WholeCell_ROIs, Soma_ROIs, Perivascular_ROIs, or Processes. Filename will be the name of the .tif file from which data was analyzed.
The parent "ROI Segmentation Analysis" folder will be located in the same directory that the .tif file was loaded from.

Data by default is saved as a .mat file. The .mat file arranges the data into a row dynamic structure, with each row containing the results of a single ROI. Scrolling across columns will reveal discrete parameter values, and nested structures of function outputs. For more information on data consolidation, tabulation, and interpretation see [Working With Data Output](#Working-With-Data-Output)


Inside the save function is a preformatted excel sheet save function for limited export of data. Should the user enable it, the resulting excel sheet contains a limited view of processed data including time series, event index, and signal properties (rise, decay, and amplitude) for the current segmentation.

### Motion Corrected Channels
If your .tif files were motion corrected during preprocessing, you may opt to export the motion corrected tifs. If selected, each corrected channel will be saved as a separate .tif containing the filename with the appended channel identification (Ex/ Chn1,Chn2, etc) inside a folder named "MotionCorrectedStacks" also located in the same directory as the raw .tif file.

Motion corrected channels will be the exported array for reduced Z series exports (See Below)

.tif
### Frame Reduction .Tif Export for Video Production
For high frame rate acquistions, videos produced from .tifs may be cumbersome in data size or exceed the frame rate capabilities of common display hardware, especially when played at accelerated speed. To remedy this a save feature was added to facilitate the export of video ready stacks with data and frame reduction. 

In the save window, checking the box labeled "Create Reduced Z Series" will open an input dialog box for setting the intended playback speed (ex/ 1x, 5x, ..., nx speed), the intended video frame rate, and data averaging factor. The data reduction rate is determined by the ratio of the product of the original data acquistion rate with desired playback speed of the file, and the target monitor frame rate.

$$ R = \frac{Acquistion Rate  *  Playback Speed}{Monitor Frame Rate} $$

Therefore, for perceptually minimal loss of individual frame information, the third input ("Temporal Resultion smoothing") value should be above the reduction ratio.


.gif options
Folders


## MATLAB Workspace Integration


### Working With .mat Data Output
Within MATLAB, opening the data .mat file will load a file's segmentation data structure into the base workspace.










## Troubleshooting and help



##
