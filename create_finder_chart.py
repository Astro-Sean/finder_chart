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

from matplotlib.patches import Polygon
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp




# =============================================================================
# Configure your Finding Chart
# =============================================================================

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

# =============================================================================
# 
# =============================================================================

def format_ra_dec(ra_deg, dec_deg):
    # Create a SkyCoord object with RA and Dec in degrees
    coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)

    # Format RA and Dec in the desired string format
    ra_str = coord.ra.to_string(unit=u.hour, sep=':', precision=3, pad=True)
    dec_str = coord.dec.to_string(sep=':', precision=2, alwayssign=True, pad=True)

    return f"{ra_str}, {dec_str}"

# =============================================================================
# 
# =============================================================================

if (title_name is None): title_name = tname


comp_image = fits.open(fpath)
image = comp_image[0].data
header = comp_image[0].header
wcs_out = WCS(header)


w, footprint = find_optimal_celestial_wcs(comp_image,auto_rotate=False)
image, _ = reproject_interp((image, header), w, shape_out=footprint)

scale_arcsec_per_pixel  = proj_plane_pixel_scales(w)[0]*3600

# =============================================================================
# 
# =============================================================================

coords = SkyCoord(ra,dec, unit = (u.deg,u.deg))
pixels = w.world_to_pixel(coords)

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

Figure will be saved to '{saveloc}' with the filename '{fname}'

"""

print(s)

# =============================================================================
# MAKE THE FINDER CHART
# =============================================================================

fig = plt.figure(figsize = (7,7))

# =============================================================================
# 
# =============================================================================

ra_offset_arcsec /=3600
dec_offset_arcsec /=3600

slit_length_arcsec/=3600
slit_width_arcsec/=3600


cutsize = cutsize/ 60 # change to degrees
window = (cutsize*3600)/scale_arcsec_per_pixel

inset_cutsize = inset_cutsize*u.arcsec
inset_size = inset_cutsize.value/scale_arcsec_per_pixel

# =============================================================================
# 
# =============================================================================

image_ax1 = Cutout2D(image.copy(),(x_pix,y_pix),(window,window),wcs = w)

w = image_ax1.wcs

image_ax2_tmp = Cutout2D(image.copy(),(x_pix,y_pix),(inset_size,inset_size),wcs = w)
image_ax2 = image_ax1



ax1 = fig.add_subplot(111,projection=image_ax1.wcs)

left, bottom, width, height = [0.66, 0.66, 0.25, 0.25]
ax2 = fig.add_axes([left, bottom, width, height],projection=image_ax2.wcs)

ax2.set_aspect('equal')

# vmin_p,vmax_p = (ZScaleInterval(nsamples = 1200)).get_limits(image_ax1.data[np.isfinite(image_ax1.data)])
vmin_p,vmax_p = np.nanpercentile(image_ax1.data[np.isfinite(image_ax1.data)], [vmin,vmax])
ax1.imshow(image_ax1.data,
           cmap='gray_r',
           origin='lower', 
           interpolation = 'None',
           vmin=vmin_p, 
           vmax=vmax_p)

vmin_p,vmax_p = (ZScaleInterval(nsamples = 1200)).get_limits(image_ax2_tmp.data[np.isfinite(image_ax2_tmp.data)])
# vmin_p,vmax_p = np.nanpercentile(image_ax2.data[np.isfinite(image_ax2.data)], [vmin,vmax])
ax2.imshow(image_ax2.data,
           cmap='gray_r',
           origin='lower',
           interpolation = 'None',
           vmin=vmin_p, 
           vmax=vmax_p)

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

target_location = format_ra_dec(ra, dec)

ra_slit = ra - ra_offset_arcsec
dec_slit = dec - dec_offset_arcsec
if ra_offset_arcsec!=0 or dec_offset_arcsec!=0:
    slit_location = format_ra_dec(ra_slit, dec_slit)

    print('\n\nTarget location ',target_location)
    print('Slit location: ',slit_location)
    
    
    
s = f"Target\n{target_location}"
if plot_slit:

    sky_coord = SkyCoord(ra=ra, dec=dec, unit='deg')
    
    # Calculate the paragalactic angle for the RA/Dec position
    sky_coord_slit = SkyCoord(ra=ra, dec=dec, unit='deg')
    gal_coord = sky_coord_slit.galactic
    paragalactic_angle = np.arctan2(gal_coord.b.degree, gal_coord.l.degree)
    
    
    
    # Calculate the paragalactic angle for the RA/Dec position
    sky_coord_slit = SkyCoord(ra=ra_slit, dec=dec_slit, unit='deg',)
    
    
    # Calculate the paragalactic angle for the RA/Dec position
    # sky_coord_slit = SkyCoord(ra=ra, dec=dec, unit='deg')
    gal_coord = sky_coord_slit.galactic
    paragalactic_angle = np.arctan2(gal_coord.b.degree, gal_coord.l.degree)
    
    paragalactic_angle = np.deg2rad(slit_angle)
    
    # Define the half-lengths for easy calculation of corners
    half_length = slit_length_arcsec / 2
    half_width = slit_width_arcsec / 2
    
    # Define corner positions based on the central point and paragalactic angle
    # Offset directions: +paragalactic_angle (along length), ±π/2 (perpendicular for width)
    corner1 = sky_coord_slit.directional_offset_by(paragalactic_angle, half_length * u.deg).directional_offset_by(paragalactic_angle + np.pi / 2, half_width * u.deg)
    corner2 = sky_coord_slit.directional_offset_by(paragalactic_angle, half_length * u.deg).directional_offset_by(paragalactic_angle - np.pi / 2, half_width * u.deg)
    corner3 = sky_coord_slit.directional_offset_by(paragalactic_angle + np.pi, half_length * u.deg).directional_offset_by(paragalactic_angle - np.pi / 2, half_width * u.deg)
    corner4 = sky_coord_slit.directional_offset_by(paragalactic_angle + np.pi, half_length * u.deg).directional_offset_by(paragalactic_angle + np.pi / 2, half_width * u.deg)
    
    # Convert corners to pixel coordinates for plotting
    corner1_pix = w.world_to_pixel(corner1)
    corner2_pix = w.world_to_pixel(corner2)
    corner3_pix = w.world_to_pixel(corner3)
    corner4_pix = w.world_to_pixel(corner4)
    
    slit_angle_str = f'{slit_angle}$\degree$' if slit_angle!=0 else 'paragalactic' # degs
    
    for ax in[ax1,ax2]:
        slit_polygon = Polygon([corner1_pix, corner2_pix, corner3_pix, corner4_pix], closed=True, edgecolor='red', facecolor='none', 
                               linewidth=0.5, label=f"{slit_length_arcsec*3600:.0f}'' $\\times$ {slit_width_arcsec*3600:.1f}'' slit\nPA: {slit_angle_str}")
        ax.add_patch(slit_polygon)

    if ra_offset_arcsec!=0 or dec_offset_arcsec!=0:
        
        s += f"\n\nSlit offset from target\n$\Delta$Ra = {ra_offset_arcsec}''\n$\Delta$Dec = {dec_offset_arcsec}''"
        
    s+=f' \n\nSlit position at PA: {slit_angle_str}'

ax1.text(0.02, 0.98, s, 
      va='top', 
      ha='left', 
      transform=ax1.transAxes, 
      color='red', 
      bbox=dict(facecolor='white', edgecolor='none'))
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

notch = 5/scale_arcsec_per_pixel

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
    fontsize = 12,
    arrowprops=dict(color='red', lw = 2.5,arrowstyle="<|-",),)


ax1.annotate(
    "E", xy=(xpos_ax1,ypos_ax1), xycoords='data',
    xytext=(xpos_ax1-length,ypos_ax1),
    ha = 'right',
    va = 'center',
    color='red',
    fontsize = 12,
    arrowprops=dict(lw = 2.5,arrowstyle="<|-",fc = 'none',ec = 'red'),)

ax1.scatter(xpos_ax1,ypos_ax1,marker = 'o',color = 'red',s = 25)


xo = +0.15
yo = +0.15

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
             fontsize = 12,
             arrowprops=dict(color='red',lw = 2.5,arrowstyle="-"))


# =============================================================================
# 
# =============================================================================

if not plot_slit:
    xpos_ax2 = x_pix_ax2 - 1 * inset_scalebar/scale_arcsec_per_pixel
    ypos_ax2 = y_pix_ax2 -1 * inset_scalebar/scale_arcsec_per_pixel
    
else:
    
    x3, y3 = corner3_pix
    x4, y4 = corner4_pix
    x_mid = (x3 + x4) / 2
    y_mid = (y3 + y4) / 2

    xpos_ax2 =   x_mid - inset_scalebar/scale_arcsec_per_pixel
    ypos_ax2 = y_mid - 0.5 * inset_scalebar/scale_arcsec_per_pixel
    

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


fig.tight_layout()

# =============================================================================
# 
# =============================================================================
plt.show(block = False)



plt.savefig(os.path.join(saveloc,fname),
            bbox_inches='tight',dpi = 300)


print('\n<> Done <>')