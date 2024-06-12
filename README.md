[arxiv_link_SN2023fyq]: https://arxiv.org/abs/2401.15148


# create_finder_chart.py

This is a relatively simple code to produce a science ready finder chart to be used to transient identification.



## Usage

To produce your own finder chart - edit the following block of code in *create_finder_chart.py*. The following code was used when producing a finder chart used in the followup of [SN2023fyq](arxiv_link_SN2023fyq)
```python

# For the filename and will be used as a label when pointing out the transient
tname = 'SN2023fyq'

# The title of the figure
title_name = 'SN2023fyq || NOT/ALFOSC || r band || 2023-12-08'

# Where you file is - WCS must be correct in the files header
fpath = '/home/seanbrennan/Desktop/SN2023fyq/images/NOT_ALFOSC/ZTF22abzzvln_NOT_ALFOSC_20231208_r_SDSS_wcs.fits'

# The location you are saving it to
saveloc = '/home/seanbrennan/Desktop/finder_charts/github'


# Ra and Dec in degress of the target
ra,dec = 186.441143074, +12.6635758235

# Size of image in arcmins
cutsize  = 5

# size of scalebar in arcmins
scalebar = 1

# # size of cutout in arcseconds
inset_cutsize = 15

# size of cutout in acrseconds
inset_scalebar = 5
```

> [!CAUTION]
> This code does not solve for WCS (world corrdnate system) values and used the transients Right Acensoin and Declincation as inputs. Make sure you have accurate WCS values written to your fits files otherwise your finder chart may be incorrect!


## Output

An example of the output of *create_finder_chart.py* is given in Figure. 1

<p align="center">
  <img src="./SN2023fyq_finder.jpeg" alt="Image" width = 600>
  <br>
  <em>Figure 1: .</em>
</p>
