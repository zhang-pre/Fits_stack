
from constants_calculate import *
from astropy.io import fits
from test import *
import numpy as np



# 设置要多少常数模型进行图像堆叠
def run(constant, filelist):
    # constant = 12

    # 设置图像的长宽
    # x_length = 2048
    # y_length = 2048
    x_length = 4112
    y_length = 4096
    # 遍历文件夹下的所有文件名
    # filelist = walkFile(filename)
    # 看文件夹里面有多少个文件
    imagemumber = filelist.__len__()
    # 以第一个文件作为主图像
    masterfile = filelist[0]

    # 定义输出的文件名及其路径
    outputpath = masterfile.split(".")[0] + "_" + str(constant) + '+%d.fit' %(imagemumber-1)
    outputpath = "result_star/" + outputpath.split("/")[1]


    x_position_master, y_position_master, mag_master, header1, data1 = culculate_master_starposition(masterfile)
    # i = 1-9
    # 再一次打开并找到主图像的头文件和数组信息
    hdu1 = fits.open(masterfile)
    data1 = hdu1[0].data
    data1 = np.squeeze(data1,axis=0)
    print("shape:")
    print(data1.shape)
    header1 = hdu1[0].header

    # 找到主图像头文件中的曝光时间
    # master_exposuretime = header1["EXPTIME"]
    master_exposuretime = header1["EXPOSURE"]
    # hdu2 = fits.open("original/N1512304921_1.fits")
    # data2 = hdu2[0].data
    #
    # hdu3 = fits.open("original/N1512305318_1.fits")
    # data3 = hdu3[0].data

    # 定义x和y的偏移量，可能是无效代码行
    biasconst_x = []
    biasconst_y = []

    # 定义一个数据集用来存储二维数组，以便叠加使用
    dataset = []
    dataset.append(data1)

    # const1 = [0.9987918212750526, 0.00035767239113580204, 1.9802960275532298, 0.002324571539825811, 0.9989914727611116, 14.247450233634613]
    # const2 = [1.0083044087349118, -0.002658483129478224, -14.40817170378288, 0.0014984565475660934, 0.9989749905640048, 30.982135861673736]

    # 这个for循环用于堆叠除主图像外其他图像的堆叠
    for i in range(1, imagemumber):
        print(i)
        # 单幅图像的文件名
        stackfile = filelist[i]
        # 计算常数模型，
        if constant == 6:
            consts, data2, eposuretime2 = culculate_six_constant(x_position_master, y_position_master, mag_master, stackfile, data1)
            # 常数模型的转换
            data2_trans = merge_plate_6model(data2, consts, x_length, y_length)
        if constant == 12:
            consts, data2, eposuretime2 = culculate_twelve_constant(x_position_master, y_position_master, mag_master, stackfile, data1)
            # 常数模型的转换
            data2_trans = merge_plate_12model(data2, consts, x_length, y_length)
        if constant == 20:
            consts, data2, eposuretime2 = culculate_twenty_constant(x_position_master, y_position_master, mag_master, stackfile, data1)
            # 常数模型的转换
            data2_trans = merge_plate_20model(data2, consts, x_length, y_length)
        if constant == 30:
            consts, data2, eposuretime2 = culculate_thirty_constant(x_position_master, y_position_master, mag_master, stackfile, data1)
            # 常数模型的转换
            data2_trans = merge_plate_30model(data2, consts, x_length, y_length)

        # 增加主图像的曝光时间
        master_exposuretime = master_exposuretime + eposuretime2
        # 数据集的叠加
        dataset.append(data2_trans)

    # data2_trans = merge_plate(data2, const1)
    # dataset.append(data2_trans)
    #

    data = np.stack(dataset, axis=0)
    data_stack = np.median(data, axis=0)
    # data_stack = data.median(axis=0, dtype="float64")
    header1["EXPTIME"] = master_exposuretime
    header1["EXPOSURE"] = master_exposuretime
    fits.writeto(outputpath, data_stack, header=header1, overwrite=True)

a = walkFile(filename)
# a=[1,2,3,4,5,6,7,8,9,10]
for i in range(0, len(a), 10):
    grouplist = a[i:i+10]
    print(grouplist)

    for j in [6, 12, 20]:
        if j==20:
            continue
        print(f"fist:{i}")
        run(j, grouplist)


