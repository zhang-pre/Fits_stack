
import numpy as np
from astropy.io import fits
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.io import fits
from scipy.interpolate import griddata
from scipy.optimize import leastsq

from test import *

# masterfile = ''
# masterfile = "original//A_99942-Ic-20s-20210414_212007_flat.fit"
# stackfile = "original//A_99942-Ic-20s-20210414_212037_flat.fit"

fwhm = 5.5
# 阈值必须设置为整数
threshold = 4
# 如果是空间可以设得比较小
deltaMag = 1.5

def culculate_master_starposition(masterfile):
    """传入主图像的文件名，输出主图像对应搜出每颗星的X,Y,粗略的仪器星等，头文件，图像二维数组
    :param masterfile:主图像的文件名
    :return:每颗星的X(列表),Y（列表）,粗略的仪器星等（列表），头文件，图像二维数组
    """
    hdu1 = fits.open(masterfile)
    data1 = hdu1[0].data
    data1 = np.squeeze(data1, axis=0)
    header1 = hdu1[0].header

    mean, median, std = sigma_clipped_stats(data1, sigma=3.0)
    background = 3 * median - 2 * mean
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(data1-background)

    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    sourceslist = list(sources)

    x_position_master = []
    y_position_master = []
    mag_master = []

    for i in range(sourceslist.__len__()):
        sourceslist[i] = str(sourceslist[i])
        sourceslist[i] = sourceslist[i].split("\n")
        sourceslist[i][2] = str(sourceslist[i][2]).split()
        x_position = float(sourceslist[i][2][1])
        x_position = round(x_position, 3)
        y_position = float(sourceslist[i][2][2])
        y_position = round(y_position, 3)
        mag = float(sourceslist[i][2][10])
        mag = round(mag, 3)
        x_position_master.append(x_position)
        y_position_master.append(y_position)
        mag_master.append(mag)

    # ---------------------------------------
    # x_position_master = [612.842, 611.822, 615.103]
    # y_position_master = [409.449, 413.612, 416.961]
    # mag_master = [-10.0, -10.0, -10.0]

    return x_position_master, y_position_master, mag_master, header1, data1

# x_position_master, y_position_master, mag_master, header1, data1 = culculate_master_starposition(masterfile)

def culculate_six_constant(x_position_master, y_position_master, mag_master, stackfile, data1):


    hdu2 = fits.open(stackfile)
    data2 = hdu2[0].data
    data2 = np.squeeze(data2,axis=0)
    header2 = hdu2[0].header
    # exposuretime = header2["EXPTIME"]
    exposuretime = header2["EXPOSURE"]
    mean, median, std = sigma_clipped_stats(data2, sigma=3.0)
    background = 3 * median - 2 * mean
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)

    sources2 = daofind(data2-background)

    cat1 = daofind(data1-background)

    cat2 = daofind(data2-background)

    x1 = cat1['xcentroid']
    y1 = cat1['ycentroid']
    x2 = cat2['xcentroid']
    y2 = cat2['ycentroid']

    ncat1 = len(cat1)
    ncat2 = len(cat2)
    XX = []
    YY = []

    for i in range(ncat2):
        XX.extend((x1 - x2[i]))
        YY.extend((y1 - y2[i]))
    XX = np.array(XX)
    YY = np.array(YY)
    xhist, xbins = np.histogram(XX, range=[-200, 200], bins=401)
    yhist, ybins = np.histogram(YY, range=[-200, 200], bins=401)

    idx = np.argmax(xhist)
    xsht0 = (xbins[idx] + xbins[idx + 1]) / 2.0
    idx = np.argmax(yhist)
    ysht0 = (ybins[idx] + ybins[idx + 1]) / 2.0

    print(xsht0, ysht0)

    mask = (np.abs(XX - xsht0) < 3) & (np.abs(YY - ysht0) < 3)
    print(mask.sum())
    xsht1 = np.median(XX[mask])
    ysht1 = np.median(YY[mask])

    print(xsht1, ysht1)

    for col in sources2.colnames:
        sources2[col].info.format = '%.8g'  # for consistent table output
    sourceslist = list(sources2)

    x_position_stack = []
    y_position_stack = []
    mag_stack = []

    for i in range(sourceslist.__len__()):
        sourceslist[i] = str(sourceslist[i])
        sourceslist[i] = sourceslist[i].split("\n")
        sourceslist[i][2] = str(sourceslist[i][2]).split()
        x_position = float(sourceslist[i][2][1])
        x_position = round(x_position, 3)
        y_position = float(sourceslist[i][2][2])
        y_position = round(y_position, 3)
        mag = float(sourceslist[i][2][10])
        mag = round(mag, 3)
        x_position_stack.append(x_position)
        y_position_stack.append(y_position)
        mag_stack.append(mag)

    x_list_1 = []
    y_list_1 = []


    x_list_n = []
    y_list_n = []

    for i in range(len(x_position_master)):
        x_1 = x_position_master[i]
        y_1 = y_position_master[i]
        mag_1 = mag_master[i]
        if mag_1 < -0.3:
            for j in range(len(x_position_stack)):
                x_2 = x_position_stack[j]
                y_2 = y_position_stack[j]
                mag_2 = mag_stack[j]

                if abs(x_1-x_2 - xsht1) < 3.0 and abs(y_1 - y_2 - ysht1) < 3.0 and abs(mag_1-mag_2) < deltaMag:
                    x_list_1 = np.append(x_list_1, [x_1])
                    y_list_1 = np.append(y_list_1, [y_1])
                    x_list_n = np.append(x_list_n, [x_2])
                    y_list_n = np.append(y_list_n, [y_2])

    print(x_list_1)
    print(x_list_n)

    def func_ksi(p1, x1, y1):
        a1, b1, c1 = p1
        return a1 * x1 + b1 * y1 + c1

    def error_ksi(p1, x1, y1, z1):
        return func_ksi(p1, x1, y1) - z1

    def func_ita(p2, x2, y2):
        a2, b2, c2 = p2
        return a2 * x2 + b2 * y2 + c2

    def error_ita(p2, x2, y2, z2):
        return func_ita(p2, x2, y2) - z2

    # 其他图像转换到第一幅图像
    p0 = [0, 0, 0]
    Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_n))
    a, b, c = Para_ksi[0]
    Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_n))
    e, f, g = Para_ita[0]
    const_six = [a, b, c, e, f, g]

    print(const_six)
    return const_six, data2, exposuretime

#有自己加的注释
# 找到恒星点，算两幅图像中恒星点的位置，计算常数模型，后面所有的数据就可以拿这个常数模型进行转换
def culculate_twelve_constant(x_position_master, y_position_master, mag_master, stackfile, data1):

    # data1是主文件，data2是准备堆叠的文件
    hdu2 = fits.open(stackfile)
    data2 = hdu2[0].data
    data2 = np.squeeze(data2,axis=0)
    header2 = hdu2[0].header
    exposuretime = header2["EXPOSURE"]
    # exposuretime = header2["EXPTIME"]
    mean, median, std = sigma_clipped_stats(data2, sigma=3.0)
    background = 3 * median - 2 * mean
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)

    sources2 = daofind(data2-background)

    cat1 = daofind(data1-background)

    cat2 = daofind(data2-background)

    x1 = cat1['xcentroid']
    y1 = cat1['ycentroid']
    x2 = cat2['xcentroid']
    y2 = cat2['ycentroid']

    ncat1 = len(cat1)
    ncat2 = len(cat2)
    XX = []
    YY = []
# 用来计算两组坐标之间的差值并将这些差值添加到XX和YY两个列表中，因为x1和x2之间下标对应的数据不一定是同一颗星的数据，所以要用x1 - x2[i]
    for i in range(ncat2):
        XX.extend((x1 - x2[i]))
        YY.extend((y1 - y2[i]))
    XX = np.array(XX)
    YY = np.array(YY)
    # 具体来说，它将XX/YY中的数值分成401个区间（bins = 401），并计算每个区间中数值的频数。
    # range = [-200, 200]指定了直方图的范围，只包括XX/YY中在 - 200到200之间的数值。
    # 结果被存储在xhist中，而区间的边界被存储在xbins中。
    xhist, xbins = np.histogram(XX, range=[-200, 200], bins=401)
    yhist, ybins = np.histogram(YY, range=[-200, 200], bins=401)

    # 找到X和Y方向的直方图中的峰值，然后计算峰值所在区间的中点，以得到X和Y方向的位移（xsht0和ysht0）
    idx = np.argmax(xhist)
    xsht0 = (xbins[idx] + xbins[idx + 1]) / 2.0
    idx = np.argmax(yhist)
    ysht0 = (ybins[idx] + ybins[idx + 1]) / 2.0

    print(xsht0, ysht0)

    mask = (np.abs(XX - xsht0) < 3) & (np.abs(YY - ysht0) < 3)
    print(mask.sum())
    # mask可能全为false啊
    xsht1 = np.median(XX[mask])
    ysht1 = np.median(YY[mask])

    print(xsht1, ysht1)

    for col in sources2.colnames:
        sources2[col].info.format = '%.8g'  # for consistent table output
    sourceslist = list(sources2)

    x_position_stack = []
    y_position_stack = []
    mag_stack = []

    for i in range(sourceslist.__len__()):
        sourceslist[i] = str(sourceslist[i])
        sourceslist[i] = sourceslist[i].split("\n")
        sourceslist[i][2] = str(sourceslist[i][2]).split()
        x_position = float(sourceslist[i][2][1])
        x_position = round(x_position, 3)
        y_position = float(sourceslist[i][2][2])
        y_position = round(y_position, 3)
        mag = float(sourceslist[i][2][10])
        mag = round(mag, 3)
        x_position_stack.append(x_position)
        y_position_stack.append(y_position)
        mag_stack.append(mag)
    # x_list_1,y_list_1 存储主图像的坐标
    x_list_1 = []
    y_list_1 = []


    x_list_n = []
    y_list_n = []
    # 主要目的是匹配两幅图像中亮度低于-0.3的星星目标，并将它们的坐标存储在不同的列表中，以便后续的坐标转换和对齐操作。
    for i in range(len(x_position_master)):
        x_1 = x_position_master[i]
        y_1 = y_position_master[i]
        mag_1 = mag_master[i]
        # 为什么是-0，3
        if mag_1 < -0.3:
            for j in range(len(x_position_stack)):
                x_2 = x_position_stack[j]
                y_2 = y_position_stack[j]
                mag_2 = mag_stack[j]

                if abs(x_1-x_2 - xsht1) < 3.0 and abs(y_1 - y_2 - ysht1) < 3.0 and abs(mag_1-mag_2) < deltaMag:
                    x_list_1 = np.append(x_list_1, [x_1])
                    y_list_1 = np.append(y_list_1, [y_1])
                    x_list_n = np.append(x_list_n, [x_2])
                    y_list_n = np.append(y_list_n, [y_2])

    print(x_list_1)
    print(x_list_n)


    def func_ksi(p1, x1, y1):
        a1, a2, a3, a4, a5, a6 = p1
        return a1 * x1 * x1 + a2 * y1 * y1 + a3 * x1 * y1 + a4 * x1 + a5 * y1 + a6

    def error_ksi(p1, x1, y1, z1):
        return func_ksi(p1, x1, y1) - z1

    def func_ita(p2, x2, y2):
        b1, b2, b3, b4, b5, b6 = p2
        return b1 * x2 * x2 + b2 * y2 * y2 + b3 * x2 * y2 + b4 * x2 + b5 * y2 + b6

    def error_ita(p2, x2, y2, z2):
        return func_ita(p2, x2, y2) - z2


    # 其他图像转换到第一幅图像
    p0 = [0, 0, 0, 0, 0, 0]
    Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_n))
    a1, a2, a3, a4, a5, a6 = Para_ksi[0]
    Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_n))
    b1, b2, b3, b4, b5, b6 = Para_ita[0]
    const_twelve = [a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6]

    print(const_twelve)
    return const_twelve, data2, exposuretime

def culculate_twenty_constant(x_position_master, y_position_master, mag_master, stackfile, data1):


    hdu2 = fits.open(stackfile)
    data2 = hdu2[0].data
    data2 = np.squeeze(data2, axis=0)
    header2 = hdu2[0].header
    # exposuretime = header2["EXPTIME"]
    exposuretime = header2["EXPOSURE"]
    mean, median, std = sigma_clipped_stats(data2, sigma=3.0)
    background = 3 * median - 2 * mean
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)

    sources2 = daofind(data2-background)

    cat1 = daofind(data1-background)

    cat2 = daofind(data2-background)

    x1 = cat1['xcentroid']
    y1 = cat1['ycentroid']
    x2 = cat2['xcentroid']
    y2 = cat2['ycentroid']

    ncat1 = len(cat1)
    ncat2 = len(cat2)
    XX = []
    YY = []

    for i in range(ncat2):
        XX.extend((x1 - x2[i]))
        YY.extend((y1 - y2[i]))
    XX = np.array(XX)
    YY = np.array(YY)
    xhist, xbins = np.histogram(XX, range=[-200, 200], bins=401)
    yhist, ybins = np.histogram(YY, range=[-200, 200], bins=401)

    idx = np.argmax(xhist)
    xsht0 = (xbins[idx] + xbins[idx + 1]) / 2.0
    idx = np.argmax(yhist)
    ysht0 = (ybins[idx] + ybins[idx + 1]) / 2.0

    print(xsht0, ysht0)

    mask = (np.abs(XX - xsht0) < 3) & (np.abs(YY - ysht0) < 3)
    print(mask.sum())
    xsht1 = np.median(XX[mask])
    ysht1 = np.median(YY[mask])

    print(xsht1, ysht1)

    for col in sources2.colnames:
        sources2[col].info.format = '%.8g'  # for consistent table output
    sourceslist = list(sources2)

    x_position_stack = []
    y_position_stack = []
    mag_stack = []

    for i in range(sourceslist.__len__()):
        sourceslist[i] = str(sourceslist[i])
        sourceslist[i] = sourceslist[i].split("\n")
        sourceslist[i][2] = str(sourceslist[i][2]).split()
        x_position = float(sourceslist[i][2][1])
        x_position = round(x_position, 3)
        y_position = float(sourceslist[i][2][2])
        y_position = round(y_position, 3)
        mag = float(sourceslist[i][2][10])
        mag = round(mag, 3)
        x_position_stack.append(x_position)
        y_position_stack.append(y_position)
        mag_stack.append(mag)

    x_list_1 = []
    y_list_1 = []


    x_list_n = []
    y_list_n = []

    for i in range(len(x_position_master)):
        x_1 = x_position_master[i]
        y_1 = y_position_master[i]
        mag_1 = mag_master[i]
        if mag_1 < -0.3:
            for j in range(len(x_position_stack)):
                x_2 = x_position_stack[j]
                y_2 = y_position_stack[j]
                mag_2 = mag_stack[j]

                if abs(x_1-x_2 - xsht1) < 3.0 and abs(y_1 - y_2 - ysht1) < 3.0 and abs(mag_1-mag_2) < deltaMag:
                    x_list_1 = np.append(x_list_1, [x_1])
                    y_list_1 = np.append(y_list_1, [y_1])
                    x_list_n = np.append(x_list_n, [x_2])
                    y_list_n = np.append(y_list_n, [y_2])

    print(x_list_1)
    print(x_list_n)


    def func_ksi(p1, x1, y1):
        a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = p1
        return a1 * x1 * x1 * x1 + a2 * y1 * y1 * y1 + a3 * x1 * x1 * y1 + a4 * x1 * y1 * y1 + \
               a5 * x1 * x1 + a6 * y1 * y1 + a7 * x1 * y1 + \
               a8 * x1 + a9 * y1 + a10

    def error_ksi(p1, x1, y1, z1):
        return func_ksi(p1, x1, y1) - z1

    def func_ita(p2, x2, y2):
        b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = p2
        return b1 * x2 * x2 * x2 + b2 * y2 * y2 * y2 + b3 * x2 * x2 * y2 + b4 * x2 * y2 * y2 + \
               b5 * x2 * x2 + b6 * y2 * y2 + b7 * x2 * y2 + \
               b8 * x2 + b9 * y2 + b10

    def error_ita(p2, x2, y2, z2):
        return func_ita(p2, x2, y2) - z2


    # 其他图像转换到第一幅图像
    p0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_n))
    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 = Para_ksi[0]
    Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_n))
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10 = Para_ita[0]
    const_twenty = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10]

    print(const_twenty)
    return const_twenty, data2, exposuretime



def culculate_thirty_constant(x_position_master, y_position_master, mag_master, stackfile, data1):


    hdu2 = fits.open(stackfile)
    data2 = hdu2[0].data
    data2 = np.squeeze(data2, axis=0)
    header2 = hdu2[0].header
    # exposuretime = header2["EXPTIME"]
    exposuretime = header2["EXPOSURE"]
    mean, median, std = sigma_clipped_stats(data2, sigma=3.0)
    background = 3 * median - 2 * mean
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)

    sources2 = daofind(data2-background)

    cat1 = daofind(data1-background)

    cat2 = daofind(data2-background)

    x1 = cat1['xcentroid']
    y1 = cat1['ycentroid']
    x2 = cat2['xcentroid']
    y2 = cat2['ycentroid']

    ncat1 = len(cat1)
    ncat2 = len(cat2)
    XX = []
    YY = []

    for i in range(ncat2):
        XX.extend((x1 - x2[i]))
        YY.extend((y1 - y2[i]))
    XX = np.array(XX)
    YY = np.array(YY)
    xhist, xbins = np.histogram(XX, range=[-200, 200], bins=401)
    yhist, ybins = np.histogram(YY, range=[-200, 200], bins=401)

    idx = np.argmax(xhist)
    xsht0 = (xbins[idx] + xbins[idx + 1]) / 2.0
    idx = np.argmax(yhist)
    ysht0 = (ybins[idx] + ybins[idx + 1]) / 2.0

    print(xsht0, ysht0)

    mask = (np.abs(XX - xsht0) < 3) & (np.abs(YY - ysht0) < 3)
    print(mask.sum())
    xsht1 = np.median(XX[mask])
    ysht1 = np.median(YY[mask])

    print(xsht1, ysht1)

    for col in sources2.colnames:
        sources2[col].info.format = '%.8g'  # for consistent table output
    sourceslist = list(sources2)

    x_position_stack = []
    y_position_stack = []
    mag_stack = []

    for i in range(sourceslist.__len__()):
        sourceslist[i] = str(sourceslist[i])
        sourceslist[i] = sourceslist[i].split("\n")
        sourceslist[i][2] = str(sourceslist[i][2]).split()
        x_position = float(sourceslist[i][2][1])
        x_position = round(x_position, 3)
        y_position = float(sourceslist[i][2][2])
        y_position = round(y_position, 3)
        mag = float(sourceslist[i][2][10])
        mag = round(mag, 3)
        x_position_stack.append(x_position)
        y_position_stack.append(y_position)
        mag_stack.append(mag)

    x_list_1 = []
    y_list_1 = []


    x_list_n = []
    y_list_n = []

    for i in range(len(x_position_master)):
        x_1 = x_position_master[i]
        y_1 = y_position_master[i]
        mag_1 = mag_master[i]
        if mag_1 < -0.3:
            for j in range(len(x_position_stack)):
                x_2 = x_position_stack[j]
                y_2 = y_position_stack[j]
                mag_2 = mag_stack[j]

                if abs(x_1-x_2 - xsht1) < 3.0 and abs(y_1 - y_2 - ysht1) < 3.0 and abs(mag_1-mag_2) < deltaMag:
                    x_list_1 = np.append(x_list_1, [x_1])
                    y_list_1 = np.append(y_list_1, [y_1])
                    x_list_n = np.append(x_list_n, [x_2])
                    y_list_n = np.append(y_list_n, [y_2])

    print(x_list_1)
    print(x_list_n)


    def func_ksi(p1, x1, y1):
        a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15 = p1
        return a1 * x1 * x1 * x1 * x1 + a2 * y1 * y1 * y1 * y1 + a3 * x1 * x1 * x1 * y1 + a4 * x1 * y1 * y1 * y1 + a5 * x1 * x1 * y1 * y1 + \
               a6 * x1 * x1 * x1 + a7 * y1 * y1 * y1 + a8 * x1 * x1 * y1 + a9 * x1 * y1 * y1 + \
               a10 * x1 * x1 + a11 * y1 * y1 + a12 * x1 * y1 + \
               a13 * x1 + a14 * y1 + a15

    def error_ksi(p1, x1, y1, z1):
        return func_ksi(p1, x1, y1) - z1

    def func_ita(p2, x2, y2):
        b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15 = p2
        return b1 * x2 * x2 * x2 * x2 + b2 * y2 * y2 * y2 * y2 + b3 * x2 * x2 * x2 * y2 + b4 * x2 * y2 * y2 * y2 + b5 * x2 * x2 * y2 * y2 + \
               b6 * x2 * x2 * x2 + b7 * y2 * y2 * y2 + b8 * x2 * x2 * y2 + b9 * x2 * y2 * y2 + \
               b10 * x2 * x2 + b11 * y2 * y2 + b12 * x2 * y2 + \
               b13 * x2 + b14 * y2 + b15

    def error_ita(p2, x2, y2, z2):
        return func_ita(p2, x2, y2) - z2


    # 其他图像转换到第一幅图像
    p0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_n))
    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15 = Para_ksi[0]
    Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_n))
    b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15 = Para_ita[0]
    const_thirty = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15]

    print(const_thirty)
    return const_thirty, data2, exposuretime