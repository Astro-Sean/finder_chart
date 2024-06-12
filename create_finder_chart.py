#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 14:55:48 2023

@author: seanbrennan
"""
import os

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from astropy.visualization import  ZScaleInterval


from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u



from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
    
# =============================================================================
# Configure your Finding Chart
# =============================================================================

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

# =============================================================================
# 
# =============================================================================

if (title_name is None): title_name = tname


comp_image = fits.open(fpath)
image = comp_image[0].data
header = comp_image[0].header
wcs_out = WCS(header)


optimal_wcs, footprint = find_optimal_celestial_wcs(comp_image,auto_rotate=False)
image, _ = reproject_interp((image, header), optimal_wcs, shape_out=footprint)

scale_arcsec_per_pixel  = proj_plane_pixel_scales(optimal_wcs)[0]*3600

# =============================================================================
# 
# =============================================================================

coords = SkyCoord(ra,dec, unit = (u.deg,u.deg))
pixels = optimal_wcs.world_to_pixel(coords)


x_pix = pixels[0]
y_pix = pixels[1]



fname = fr'{tname}_finder.jpeg'

# =============================================================================
# 
# =============================================================================

s = f""" \n\n> Creating a finding chart for {tname} < 


Source of interest located at x = {x_pix:.1f}, y = {y_pix:.1f}

Image will be {cutsize:.1f} arcmins wide wtih a {inset_cutsize:.1f} arcsecs wide inset

Figure will have the title '{title_name}' 

Figure will be saved to {saveloc} with the filename '{fname}'

"""

print(s)

# =============================================================================
# MAKE THE FINDER CHART
# =============================================================================

fig = plt.figure(figsize = (7,7))

# =============================================================================
# 
# =============================================================================

cutsize = cutsize/ 60 # change to degrees
window = (cutsize*3600)/scale_arcsec_per_pixel

inset_cutsize = inset_cutsize*u.arcsec
inset_size = inset_cutsize.value/scale_arcsec_per_pixel

# =============================================================================
# 
# =============================================================================

image_ax1 = Cutout2D(image.copy(),(x_pix,y_pix),(window,window),wcs = optimal_wcs)


image_ax2_tmp = Cutout2D(image.copy(),(x_pix,y_pix),(inset_size,inset_size),wcs = optimal_wcs)
image_ax2 = image_ax1



ax1 = fig.add_subplot(111,projection=image_ax1.wcs)

left, bottom, width, height = [0.65, 0.65, 0.25, 0.25]
ax2 = fig.add_axes([left, bottom, width, height],projection=image_ax2.wcs)


vmin,vmax = (ZScaleInterval(nsamples = 1200)).get_limits(image_ax1.data[np.isfinite(image_ax1.data)])
vmin,vmax = np.nanpercentile(image_ax1.data[np.isfinite(image_ax1.data)], [0.5,99.9])
ax1.imshow(image_ax1.data,
           cmap='gray_r',
           origin='lower', 
           vmin=vmin, 
           vmax=vmax)

vmin,vmax = (ZScaleInterval(nsamples = 1200)).get_limits(image_ax2_tmp.data[np.isfinite(image_ax2_tmp.data)])
ax2.imshow(image_ax2.data,
           cmap='gray_r',
           origin='lower',
           vmin=vmin, 
           vmax=vmax)


lon = ax2.coords[0]   
lat = ax2.coords[1]

lat.set_ticklabel_visible(False)
lon.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
lon.set_ticks_visible(False)


lon = ax1.coords[0]
lat = ax1.coords[1]

lon.set_axislabel('Right Ascension' )
lat.set_axislabel('Declination')

ax1.set_title(title_name)

# =============================================================================
# 
# =============================================================================

ax2.xaxis.set_ticks([])
ax2.yaxis.set_ticks([])


# =============================================================================
# 
# =============================================================================


coords = SkyCoord(ra,dec, unit = (u.deg,u.deg))

pixels = image_ax1.wcs.world_to_pixel(coords)
x_pix_ax1 = pixels[0]
y_pix_ax1 = pixels[1]

pixels = image_ax2.wcs.world_to_pixel(coords)
x_pix_ax2 = pixels[0]
y_pix_ax2 = pixels[1]


# =============================================================================
# 
# =============================================================================


notch = 10/scale_arcsec_per_pixel


ax1.annotate('',xy=(x_pix_ax1,y_pix_ax1-notch/5), xycoords='data',
    xytext=(x_pix_ax2,y_pix_ax1-notch-notch/5),
    ha = 'center',
    va = 'bottom',
    color='red',
    arrowprops=dict(color='red', lw = 1.5,arrowstyle="-",),)



ax1.annotate('',xy=(x_pix_ax1-notch/5,y_pix_ax1), xycoords='data',
    xytext=(x_pix_ax1-notch/5-notch,y_pix_ax1),
    ha = 'center',
    va = 'bottom',
    color='red',
    arrowprops=dict(color='red', lw = 1.5,arrowstyle="-",),)

# =============================================================================
# 
# =============================================================================
fov_deg = inset_cutsize.to('deg').value


notch = 3/scale_arcsec_per_pixel



# =============================================================================
# 
# =============================================================================

xpos_ax1 = ax1.get_xlim()[1] - (ax1.get_xlim()[1] - ax1.get_xlim()[0])*0.1
ypos_ax1 = ax1.get_ylim()[0] + (ax1.get_ylim()[1] - ax1.get_ylim()[0])*0.1

length = (ax1.get_xlim()[1] - ax1.get_xlim()[0])*0.2

ax1.annotate(
    "N", xy=(xpos_ax1,ypos_ax1), xycoords='data',
    xytext=(xpos_ax1,ypos_ax1+length),
    ha = 'center',
    va = 'bottom',
    color='red',
    arrowprops=dict(color='red', lw = 1.5,arrowstyle="<|-",),)


ax1.annotate(
    "E", xy=(xpos_ax1,ypos_ax1), xycoords='data',
    xytext=(xpos_ax1-length,ypos_ax1),
    ha = 'right',
    va = 'center',
    color='red',
    arrowprops=dict(lw = 1.5,arrowstyle="<|-",fc = 'none',ec = 'red'),)

ax1.scatter(xpos_ax1,ypos_ax1,marker = 'o',color = 'red',s = 25)




xo = +1
yo = +1

ax2.annotate(tname,xy=(x_pix_ax2+xo,y_pix_ax2+yo), xycoords='data',
    xytext=(x_pix_ax2+notch*0.8,y_pix_ax2+notch*2),
    ha = 'center',
    va = 'bottom',
    color='red',
    arrowprops=dict(color='red', lw = 0.5,arrowstyle="-|>",),)

# =============================================================================
# 
# =============================================================================

xpos = ax1.get_xlim()[0] + (ax1.get_xlim()[1] - ax1.get_xlim()[0])*0.1
ypos = ax1.get_ylim()[0] + (ax1.get_ylim()[1] - ax1.get_ylim()[0])*0.1

def format_value(value):
    if value.is_integer():
        return int(value)
    else:
        return round(value, 1)



scalebar= scalebar* 60
if scalebar < 60:
    label = str()+"''"
elif scalebar < 360 and scalebar >= 60:
    # scale = [int(scale/60) if scale/600 ==0 else float(scale/60)]
    label = str(format_value(scalebar/60))+"'"
elif scalebar > 360:
    label = str(format_value(scalebar/360))+"d"
    
ax1.annotate(label,
             xy=(xpos,ypos), xycoords='data',
             xytext=(xpos+scalebar/scale_arcsec_per_pixel,ypos),
             ha = 'left',
             va = 'center',
             color='red',
             fontsize = 9,
             arrowprops=dict(color='red',lw = 1.5,arrowstyle="-"))


# =============================================================================
# 
# =============================================================================
xpos_ax2 = x_pix_ax2 - 1 * inset_scalebar/scale_arcsec_per_pixel
ypos_ax2 = y_pix_ax2 -1 * inset_scalebar/scale_arcsec_per_pixel

if inset_scalebar < 60:
    label = str(inset_scalebar)+"''"
elif inset_scalebar < 360 and inset_scalebar > 60:
    # scale = [int(scale/60) if scale/600 ==0 else float(scale/60)]
    label = str(format_value(inset_scalebar/60))+"'"
elif inset_scalebar > 360:
    label = str(format_value(inset_scalebar/360))+"d"
    
ax2.annotate(label,
             xy=(xpos_ax2,ypos_ax2), xycoords='data',
             xytext=(xpos_ax2+inset_scalebar/scale_arcsec_per_pixel,ypos_ax2),
             ha = 'left',
             va = 'center',
             color='red',
             fontsize = 6,
             arrowprops=dict(color='red',lw = 1.5,arrowstyle="-"))



# =============================================================================
# 
# =============================================================================



ax2.set(xlim = [x_pix_ax2 - inset_size-1,x_pix_ax2+inset_size-1],
        ylim = [y_pix_ax2 - inset_size-1,y_pix_ax2+inset_size-1])


ax1.indicate_inset_zoom(ax2, edgecolor="red")

# ax2.spines['bottom'].set_color('red')
# ax2.spines['top'].set_color('red')
# ax2.spines['left'].set_color('red')
# ax2.spines['right'].set_color('red')
# ax2.axis('off')



fig.tight_layout()
# =============================================================================
# 
# =============================================================================
plt.show(block = False)



plt.savefig(os.path.join(saveloc,fname),
            bbox_inches='tight',dpi = 300)


print('> Done <')