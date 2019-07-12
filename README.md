System Requirements
===================

* Compiled and tested with MATLAB R2019a, R2017b on Windows 10 64-bit

Installation Guide
==================

To download the source code:

1. Click on the Repository tab above.
2. Click on the Download icon and select the desired format. If necessary, unzip the files.
3. Add the path to the downloaded directory on MATLAB. The easiest way is to right-click on the folder in the MATLAB IDE, then select **Add to Path** > **Selected Folders and Subfolders**.

Typical install time: 1 minute

Usage
=====

The code is divided into two main functions: (1) segmentation and subsequent tracking of cells, and (2) data anlaysis.

Segmentation and tracking
-------------------------

1. Create a new FRETtracker object:

  F = FRETtracker;
  
2. Run the tracking using the default parameters:

  processFiles(F, 'example.nd2')

The code will output:
* A MAT-file (*.mat) that contains the tracked data
* An AVI-file (*.avi) that displays the segmentation and tracking results. In this movie, the cell outlines are shown in green, while the tracks are shown as white lines. The numbers next to each cell correspond to their ID.

Expected run time for tracking code: Approximately 1 - 2 hours per movie.

Data analysis
-------------

1. Load the MAT-file created by the tracking code. This should create a new variable in the Workspace called 'trackData'

  load filename.mat

2. Retrieve data from a single track

  track1 = getTrack(trackData, 1);
  
3. Retrieve the CDK2 intensities from the nucleus (nucl) and the cytoplasm (cyto):

  CDK2nucl = getData(track1, 'CDK2nucl');
  CDK2cyto = getData(track1, 'CDK2cyto');
  
4. Compute the CDK2 ratio (nucl/cyto) and plot

  plot(CDK2nucl ./ CDK2ctyo)
  
For additional examples and a fullusage guide, please see the [Wiki](https://biof-git.colorado.edu/biofrontiers-imaging/palmer-zinc-cell-cycle/wikis/home)