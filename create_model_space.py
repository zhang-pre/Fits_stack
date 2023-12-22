# from constants_calculate import *
from astropy.io import fits
# from test import *
import numpy as np
import os

filename = "model/"

def walkFile(filename):
    for root, dirs, files in os.walk(filename):

        # root 表示当前正在访问的文件夹路径
        # dirs 表示该文件夹下的子目录名list
        # files 表示该文件夹下的文件list
        filelist = []
        # 遍历文件
        for f in files:
            print(os.path.join(root, f))
            filelist.append(os.path.join(root, f))
    return filelist


imagemumber = 7
filelist = walkFile(filename)
masterfile = filelist[0]

outputpath = masterfile.split(".")[0] + '+%d_model.fit' %(imagemumber-1)
outputpath = "result_star/" + outputpath.split("/")[1]

length = filelist.__len__()
dataset = []
header = ""
for i in range(0, length):
    print(i)
    hdu1 = fits.open(filelist[i])
    data1 = hdu1[0].data
    header1 = hdu1[0].header
    header = header1
    dataset.append(data1)

data = np.stack(dataset, axis=0)
data_stack = np.median(data, axis=0)
# data_stack = data.mean(axis=0)
fits.writeto(outputpath, data_stack, header=header, overwrite=True)

outputpath2 = "model_substract/" + outputpath.split("/")[1]

data_new = fits.open("result_star/N1512304524_1+6.fit")[0].data
headernew = fits.open("result_star/N1512304524_1+6.fit")[0].header

data_substracted = data_new - data_stack + 20.0
fits.writeto(outputpath2, data_substracted, header=headernew, overwrite=True)



