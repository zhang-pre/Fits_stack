
from constants_calculate import *
from astropy.io import fits
from test import *

imagemumber = 10
filelist = walkFile(filename)
masterfile = filelist[0]

outputpath = masterfile.split(".")[0] + '+%d.fit' %(imagemumber-1)
outputpath = "result_object/" + outputpath.split("/")[1]

x_position_master, y_position_master, mag_master, header1, data1 = culculate_master_starposition(masterfile)
# i = 1-9
biasconst_x = []
biasconst_y = []
# note biasconst.txt有没有可能是天然卫星的x，y的偏移距离
f_bias = open("constbias.txt", "r+")
biaslist = f_bias.readlines()
for line in biaslist:
    line = line.split("\t")
    x = float(line[0])
    y = float(line[1])
    bias_x = biasconst_x.append(x)
    bias_y = biasconst_y.append(y)

print(biasconst_x)
print(biasconst_y)

dataset = []
dataset.append(data1)

for i in range(1, imagemumber):
    print(i)
    stackfile = filelist[i]
    consts, data2 = culculate_six_constant(x_position_master, y_position_master, mag_master, stackfile)
    print(consts)
    consts[2] = consts[2] + biasconst_x[i]
    consts[5] = consts[5] - biasconst_y[i]
    print(consts)
    data2_trans = merge_plate(data2, consts)
    dataset.append(data2_trans)

data = np.stack(dataset, axis=0)
data_stack = data.mean(axis=0, dtype="float64")

fits.writeto(outputpath, data_stack, header=header1, overwrite=True)







