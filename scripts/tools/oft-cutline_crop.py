#!/usr/bin/python

#This file is part of Open Foris Geospatial Toolkit which is free software.
#You can redistribute it and/or modify it under the terms of the 
#GNU General Public License as published by the Free Software Foundation, 
#either version 3 of the License, or (at your option) any later version.
# Visit http://www.openforis.org/tools/geospatial-toolkit.html

#Open Foris Geospatial Toolkit is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with the Open Foris Geospatial Toolkit  
#If not, see <http://www.gnu.org/licenses/>.

#####################################################################################
# Script by remi.dannunzio@fao.org
#
# Purpose : Crop a raster file to the extent of a shapefile with strict alignment 
# of the original raster and cropped output
#####################################################################################

version="1.0"
update ="17/01/2016" 

import sys
import subprocess
import os
import tempfile
import re
import shutil

try:
    from osgeo import gdal
except ImportError:
    import gdal

gdal.AllRegister()

argv = None

if argv is None:
    argv = sys.argv
    argv = gdal.GeneralCmdLineProcessor( argv )
    if argv is None:
        sys.exit( 0 )

###############################################################
##################### Test if all necessary arguments are in
###############################################################
if len(argv) != 9:
    print "\nVersion %s Last update on %s\nCrop a raster to the extent of a shapefile with alignment of original pixels\n" % (version,update)
    print "Usage: oft-cutline <-v vector_input> <-i mask_raster_input> <-o rasterized_output> [-a attribute]\n"
    sys.exit( 0 )


###############################################################
##################### Create a temp directory
###############################################################
tmpdir = tempfile.mkdtemp()
homedir = os.getcwd()


###############################################################
##################### Read the arguments passed to the command
###############################################################
i=1
while i < len(argv):
    arg = argv[i]

    if arg == '-a':
        i = i + 1
        attr = argv[i];
    
    elif arg == '-o':
        i = i + 1
        outfile = argv[i];
        
    elif arg == '-i':
        i = i + 1
        infile = argv[i];
               
    elif arg == '-v':
        i = i + 1
        vector = argv[i];
        
    elif arg[:1] == '-':
        print('Unrecognised option: %s' % arg)
        Usage()
        sys.exit( 1 )
    i=i+1

###############################################################
##################### Create basenames 
###############################################################
in_base = os.path.basename(infile)
vc_base = os.path.basename(vector)
layer   = vc_base[0:len(vc_base)-4]
outdir  = os.path.split(os.path.abspath(outfile))[0]


###############################################################
##################### Extract extent of image and pixel size
###############################################################
ul = subprocess.check_output("gdalinfo %s | grep \"Upper Left\""  % infile, shell=True)
lr = subprocess.check_output("gdalinfo %s | grep \"Lower Right\"" % infile, shell=True)
px = subprocess.check_output("gdalinfo %s | grep \"Pixel Size\""  % infile, shell=True)

im_xmin = float(ul.split("(")[1].split(",")[0])
im_ymax = float(re.sub(r'\)', '', ul.split("(")[1].split(",")[1]))
im_xmax = float(lr.split("(")[1].split(",")[0])
im_ymin = float(re.sub(r'\)', '', lr.split("(")[1].split(",")[1]))
pxsz    = float(px.split("(")[1].split(",")[0])

print "-"*40
print "Infile:  %s\nxmin : %r \nymax : %r \nxmax : %r \nymin : %r \nsize : %r\n" % (in_base,im_xmin,im_ymax,im_xmax,im_ymin,pxsz)


###############################################################
##################### Extract extent of vector
###############################################################
vc = subprocess.check_output("ogrinfo -al -so %s | grep \"Extent\"" % vector, shell=True)
vc=re.sub(r'Extent:\s','',vc)

vc_xmin = float(re.sub(r'\(', '', vc.split(") - (")[0].split(",")[0]))
vc_ymin = float(re.sub(r'\)', '', vc.split(") - (")[0].split(",")[1]))
vc_xmax = float(re.sub(r'\(', '', vc.split(") - (")[1].split(",")[0]))
vc_ymax = float(re.sub(r'\)', '', vc.split(") - (")[1].split(",")[1]))

print "-"*40
print "Mask:    %s\nvc_xmin : %r \nvc_ymax : %r \nvc_xmax : %r \nvc_ymin : %r\n" % (vc_base,vc_xmin,vc_ymax,vc_xmax,vc_ymin)


################################################################
##################### Calculate starting point for rasterization
################################################################
offset_x = int((vc_xmin-im_xmin)/pxsz)-1
offset_y = int((vc_ymin-im_ymin)/pxsz)-1

xmin = im_xmin+offset_x*pxsz
ymin = im_ymin+offset_y*pxsz


################################################################
#################### Calculate size of extent to crop
################################################################
size_x = int((vc_xmax-vc_xmin)/pxsz)+3;
size_y = int((vc_ymax-vc_ymin)/pxsz)+3;

xmax = xmin+size_x*pxsz;
ymax = ymin+size_y*pxsz;


################################################################
#################### Perform the OFGT commands
################################################################

################ Rasterize the vector with at least one pixel encompassing the border
rasterize = "gdal_rasterize -a %s -l %s -te %r %r %r %r -tr %r %r -co \"COMPRESS=LZW\" %s %s/tmp_mask.tif\n" % (attr,layer,im_xmin,im_ymin,im_xmax,im_ymax,pxsz,pxsz,vector,tmpdir)
print "-"*40 + "\n" + rasterize
os.system(rasterize)

################# Crop the original image to the same extent
translate = "gdal_translate -co \"COMPRESS=LZW\" -projwin %r %r %r %r %s  %s/tmp_crop.tif\n" % (im_xmin,im_ymax,im_xmax,im_ymin,infile,tmpdir)
print "-"*40 + "\n" + translate
os.system(translate)

################# Mask out the values from the rasterized vector mask
mask_out  = "oft-calc -um %s/tmp_mask.tif %s/tmp_crop.tif %s/tmp_out.tif <<eof \n 1 \n #1 \n eof\n" % (tmpdir,tmpdir,tmpdir)
print "-"*40 + "\n" + mask_out
os.system(mask_out)

################# Compress and output
compress = "gdal_translate -projwin %r %r %r %r -co \"COMPRESS=LZW\" %s/tmp_out.tif %s" % (vc_xmin,vc_ymax,vc_xmax,vc_ymin,tmpdir,outfile) 
print "-"*40 + "\n" + compress
os.system(compress)


################################################################
#################### Clean the temporary directory
################################################################
shutil.rmtree(tmpdir)
