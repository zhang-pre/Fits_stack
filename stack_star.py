
from constants_calculate import *
from astropy.io import fits
from test import *
import numpy as np

filelist = walkFile(filename)
imagemumber = filelist.__len__()
masterfile = filelist[0]

outputpath = masterfile.split(".")[0] + '+%d.fit' %(imagemumber-1)
outputpath = "result_star/" + outputpath.split("/")[1]

x_position_master, y_position_master, mag_master, header1, data1 = culculate_master_starposition(masterfile)
# i = 1-9

hdu1 = fits.open(masterfile)
data1 = hdu1[0].data
header1 = hdu1[0].header

print(data1)

master_exposuretime = header1["EXPTIME"]
print(master_exposuretime)
# hdu2 = fits.open("original/N1512304921_1.fits")
# data2 = hdu2[0].data
#
# hdu3 = fits.open("original/N1512305318_1.fits")
# data3 = hdu3[0].data


biasconst_x = []
biasconst_y = []

dataset = []
dataset.append(data1)

# const1 = [0.9987918212750526, 0.00035767239113580204, 1.9802960275532298, 0.002324571539825811, 0.9989914727611116, 14.247450233634613]
# const2 = [1.0083044087349118, -0.002658483129478224, -14.40817170378288, 0.0014984565475660934, 0.9989749905640048, 30.982135861673736]

for i in range(1, imagemumber):
    print(i)

    stackfile = filelist[i]


    consts, data2, eposuretime2 = culculate_six_constant(x_position_master, y_position_master, mag_master, stackfile, data1)
    master_exposuretime = master_exposuretime + eposuretime2
    data2_trans = merge_plate(data2, consts)
    dataset.append(data2_trans)

# data2_trans = merge_plate(data2, const1)
# dataset.append(data2_trans)
#

data = np.stack(dataset, axis=0)
# data_stack = np.median(data, axis=0)
data_stack = data.mean(axis=0)
header1["EXPTIME"] = master_exposuretime
header1["EXPOSURE"] = master_exposuretime
fits.writeto(outputpath, data_stack, header=header1, overwrite=True)







