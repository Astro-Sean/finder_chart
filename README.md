
# create_finder_chart.py

This is a relatively simple code using Astropy to produce a science-ready finder chart to be used for transient identification. The code will take a standard FITS image, align the image to celestial North, and annotate and highlight a given source.

This script is meant to replace the lack of a dedicated means to rapidly produce science-ready finder charts.

>  [!IMPORTANT]
> This code does not solve for WCS (World Coordinate System) values and uses the transient's Right Ascension and Declination as inputs. Make sure you have accurate WCS values written to your FITS files; otherwise, your finder chart may be incorrect!"

## Usage

To produce your own finder chart - edit the following block of code in create_finder_chart.py. The following code was used when producing a finder chart used in the follow-up of [SN 2023fyq](https://arxiv.org/abs/2401.15148).

```python
# For the filename
tname = 'SN2023fyq'

# The title of the figure
title_name = 'SN2023fyq || NOT/ALFOSC || r band || 2023-12-08'

# Where you file is
fpath = '/home/seanbrennan/Desktop/SN2023fyq/images/NOT_ALFOSC/ZTF22abzzvln_NOT_ALFOSC_20231208_r_SDSS_wcs.fits'

# The location you are saving it to
saveloc = '/home/seanbrennan/Desktop/'


# Ra and Dec in degress of the target
ra,dec = 186.441143074, +12.6635758235


# size of cutout in acrseconds
scale_cutout = 3

# Size of image in arcmins
cutsize  = 3

# size of scalebar in arcmins
scalebar = 1

# Size of inset cutout in arcseconds
inset_cutsize = 10

# Size of inset scalebar in arcseconds
inset_scalebar = 5

# lower,upper percentile values used for vmin,vmax when plotting image
# Adjust values to improve contrast
vmin,vmax = 0.5,99.5

# Plot a rectangle on the image to represent the spectrograh's slit location
plot_slit = True

# For plotting the slit - do you want the slit offset from the transient location
ra_offset_arcsec = 0
dec_offset_arcsec = 0

slit_angle = 35 # in degrees with respects to the paragalatic angle (set to 0 for paragalactic angle)

# This is just for plotting and while likely change depending on your telescope & spectrograph
slit_length_arcsec = 11  
slit_width_arcsec = 0.5   

```

Once you updated this parameters, run the code to produce your fidner chart.

```bash

 > python create_finder_chart.py

```


## Output

An example of the output of *create_finder_chart.py* is given in Figure. 1

<p align="center">
  <img src="./SN2023fyq_finder.jpeg" alt="Image" width = 600>
  <br>
  <em>Figure 1.</em>
</p>
