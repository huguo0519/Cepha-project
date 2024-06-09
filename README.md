# automatic analysis of bacteria data
# 1, read the image stack
# 2, split channels
# 3, generate a mask for channel 1 to calculate background, in a function
# 4, repeat for channel 2
# 5, before RUN, change the directory of the image folder (in line 43) and the location of csv (in line 105) and change all "\" to "/"
# 6, before RUN, change the file type in line 54

from ij import IJ, ImagePlus, ImageStack
from ij.process import ImageStatistics as IS
from ij.process import AutoThresholder
import os, csv
from loci.plugins import BF, Slicer
from loci.plugins.in import ImporterOptions
from loci.common import Region

options = IS.MEAN | IS.MEDIAN | IS.MIN_MAX

#def getStatistics(imp):
#	""" return statistics for the given imageplus """
#	IJ.run(imp, "Gaussian Blur...", "sigma=3")
#	ip = imp.getProcessor()
#	ip.setAutoThreshold("Otsu dark no-reset")
#	stats = IS.getStatistics(ip, options, imp.getCalibration())
#	return stats.mean, stats.median, stats.min, stats.max

def maskGen(imp, flag):
	""" returning the mask generated from input imp """
	""" flag = 1 means signal, flag = 0 means background"""
	IJ.run(imp, "Gaussian Blur...", "sigma=2.5")
	ip = imp.getProcessor()
	ip.setAutoThreshold("Otsu dark no-reset")
	imp.show()
	IJ.run("Convert to Mask")
	if flag == 0:
		IJ.run("Dilate")
		IJ.run("Dilate")
		IJ.run("Dilate")
		IJ.run("Dilate")
		IJ.run("Dilate")
	return imp

# folder to read all images from:
folder = "Folder Path"
# Replace Folder Path to the actual used one

# get statistics for each image in the folder
# whose file extension is ".tif"

image_stat = []
filenamelist = []
k = 0

for filename in os.listdir(folder):
	if filename.endswith(".nd2"):
		filenamelist.append(filename)
		print "Processing", filename
		file = os.path.join(folder,filename)
		
		# read in and display ImagePlus(es) with arguments
		options = ImporterOptions()
		options.setSplitChannels(1)
		options.setId(file)
		
		imps_mask_sig = BF.openImagePlus(options) # reading the file for mask generation
		imps_mask_bg = BF.openImagePlus(options) # reading the file for mask generation
		imps_imp = BF.openImagePlus(options) # reading the file again for calculation
		current_pair_total = []
		current_pair_mean = []
		for i in xrange(len(imps_mask_sig)):
			mask_sig = maskGen(imps_mask_sig[i],1)
			mask_bg = maskGen(imps_mask_bg[i],0)
			imp = imps_imp[i]
			
			ip_mask_sig = mask_sig.getProcessor().convertToFloat()
			pixels_mask_sig = ip_mask_sig.getPixels()
			ip_mask_bg = mask_bg.getProcessor().convertToFloat()
			pixels_mask_bg = ip_mask_bg.getPixels()
			ip_imp = imp.getProcessor().convertToFloat()
			pixels_imp = ip_imp.getPixels()
			above = filter(lambda i: pixels_mask_sig[i] > 0, xrange(len(pixels_mask_sig)))
			below = filter(lambda i: pixels_mask_bg[i] < 1, xrange(len(pixels_mask_bg)))
			
			background = []
			for j in xrange(len(below)):
				background.append(pixels_imp[below[j]])
			mbg = sum(background)/len(below)	# the average background value
			
			signal = []
			for j in xrange(len(above)):
				signal.append(pixels_imp[above[j]]-mbg)
				
			current_pair_total.append(sum(signal))
			current_pair_mean.append(sum(signal)/len(above))
			print "total signal:", current_pair_total[i], "mean signal:", current_pair_mean[i]
						
			k += 1
		image_stat.append([current_pair_total, current_pair_mean])

	else:
		print "Ignoring", filename
	
	
IJ.run("Close All")
print "The End of Analysis"
print image_stat
with open('Folder Path/File Name.csv', 'wb') as csvfile:  
# replace Folder Path and File Name to actual used ones 
  w = csv.writer(csvfile, delimiter=',', quotechar="\"",  
                 quoting=csv.QUOTE_NONNUMERIC)  
  w.writerow(['filename','green_total','red_total','green_mean','red_mean', 'GFP/RFP'])
  for i in xrange(len(filenamelist)):
  	current_row = image_stat[i]
  	current_total = current_row[0]
  	current_mean = current_row[1]
  	w.writerow([filenamelist[i],current_total[0], current_total[1], current_mean[0], current_mean[1], current_mean[0]/current_mean[1]])
