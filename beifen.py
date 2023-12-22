from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from astropy.stats import sigma_clipped_stats
import numpy as np
from astropy.io import fits
from scipy.interpolate import griddata
from test import *

hdu1 = fits.open("original//A_99942-Ic-20s-20210414_212007_flat.fit")
data1 = hdu1[0].data
header1 = hdu1[0].header

mean, median, std = sigma_clipped_stats(data1, sigma=3.0)
background = 3 * median - 2 * mean
daofind = DAOStarFinder(fwhm=4.0, threshold=10.0 * std)
sources = daofind(data1-background)

for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
sourceslist = list(sources)

f1 = open("positions//1-position.txt", "w", encoding='utf-8')
for i in range(sourceslist.__len__()):
    string = ''
    sourceslist[i] = str(sourceslist[i])
    sourceslist[i] = sourceslist[i].split("\n")
    sourceslist[i][2] = str(sourceslist[i][2]).split()
    x_position = float(sourceslist[i][2][1])
    x_position = round(x_position, 3)
    y_position = float(sourceslist[i][2][2])
    y_position = round(y_position, 3)
    mag = float(sourceslist[i][2][10])
    mag = round(mag, 3)
    string = str(x_position) + '   ' + str(y_position) + '   ' + str(mag) + '\n'
    f1.write(string)
f1.close()


hdu2 = fits.open("original//A_99942-Ic-20s-20210414_212037_flat.fit")
data2 = hdu2[0].data
header2 = hdu2[0].header

mean, median, std = sigma_clipped_stats(data2, sigma=3.0)
background = 3 * median - 2 * mean
daofind = DAOStarFinder(fwhm=4.0, threshold=10.0 * std)
sources = daofind(data2-background)

for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
sourceslist = list(sources)

f1 = open("positions//2-position.txt", "w", encoding='utf-8')
for i in range(sourceslist.__len__()):
    string = ''
    sourceslist[i] = str(sourceslist[i])
    sourceslist[i] = sourceslist[i].split("\n")
    sourceslist[i][2] = str(sourceslist[i][2]).split()
    x_position = float(sourceslist[i][2][1])
    x_position = round(x_position, 3)
    y_position = float(sourceslist[i][2][2])
    y_position = round(y_position, 3)
    mag = float(sourceslist[i][2][10])
    mag = round(mag, 3)
    string = str(x_position) + '   ' + str(y_position) + '   ' + str(mag) + '\n'
    f1.write(string)
f1.close()

hdu3 = fits.open("original//A_99942-Ic-20s-20210414_212107_flat.fit")
data3 = hdu3[0].data
header3 = hdu3[0].header

mean, median, std = sigma_clipped_stats(data3, sigma=3.0)
background = 3 * median - 2 * mean
daofind = DAOStarFinder(fwhm=4.0, threshold=10.0 * std)
sources = daofind(data3-background)

for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
sourceslist = list(sources)

f1 = open("positions//3-position.txt", "w", encoding='utf-8')
for i in range(sourceslist.__len__()):
    string = ''
    sourceslist[i] = str(sourceslist[i])
    sourceslist[i] = sourceslist[i].split("\n")
    sourceslist[i][2] = str(sourceslist[i][2]).split()
    x_position = float(sourceslist[i][2][1])
    x_position = round(x_position, 3)
    y_position = float(sourceslist[i][2][2])
    y_position = round(y_position, 3)
    mag = float(sourceslist[i][2][10])
    mag = round(mag, 3)
    string = str(x_position) + '   ' + str(y_position) + '   ' + str(mag) + '\n'
    f1.write(string)
f1.close()


hdu4 = fits.open("original//A_99942-Ic-20s-20210414_212157_flat.fit")
data4 = hdu4[0].data
header4 = hdu4[0].header

mean, median, std = sigma_clipped_stats(data4, sigma=3.0)
background = 3 * median - 2 * mean
daofind = DAOStarFinder(fwhm=4.0, threshold=10.0 * std)
sources = daofind(data4-background)

for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output
sourceslist = list(sources)

f1 = open("positions//4-position.txt", "w", encoding='utf-8')
for i in range(sourceslist.__len__()):
    string = ''
    sourceslist[i] = str(sourceslist[i])
    sourceslist[i] = sourceslist[i].split("\n")
    sourceslist[i][2] = str(sourceslist[i][2]).split()
    x_position = float(sourceslist[i][2][1])
    x_position = round(x_position, 3)
    y_position = float(sourceslist[i][2][2])
    y_position = round(y_position, 3)
    mag = float(sourceslist[i][2][10])
    mag = round(mag, 3)
    string = str(x_position) + '   ' + str(y_position) + '   ' + str(mag) + '\n'
    f1.write(string)
f1.close()

length = len(data2)

consts = [0.9999651436708159, -9.253577622861439e-05, -0.06014983076441392, 9.564055522149112e-05, 0.9999253414876268, 0.697723230356886]
data2_2 = merge_plate(data2, consts)
# note 线性插值
# trans = np.zeros(length * length).reshape(length, length)
#
# for i in range(length):
#     print(i)
#     for j in range(length):
#
#         x_2 = i * 1.0000 + j * 9.2542e-05 + 0.0600
#         y_2 = i * -9.5648e-05 + j * 1.0000 + (-0.6977)
#
#         x2_integer_pre = int(x_2)
#         y2_integer_pre = int(y_2)
#
#         x2_integer_post = x2_integer_pre + 1
#         y2_integer_post = y2_integer_pre + 1
#
#         if x2_integer_pre <= 2046 and x2_integer_pre >=0 and y2_integer_pre <= 2046 and y2_integer_pre >=0:
#
#             valuelist = [data2[x2_integer_pre][y2_integer_pre], data2[x2_integer_pre][y2_integer_post], data2[x2_integer_post][y2_integer_pre], data2[x2_integer_post][y2_integer_post]]
#             position = np.array([[x2_integer_pre, y2_integer_pre], [x2_integer_pre, y2_integer_post], [x2_integer_post, y2_integer_pre], [x2_integer_post, y2_integer_post]])
#             intervalue = griddata(position, valuelist, (i, j), method='linear')
#             trans[i][j] = intervalue
#
# print(trans)

data = np.stack((data1, data2_2), axis=0)

data_stack = data.mean(axis=0, dtype="float64")

print(data_stack)
# print(header1)


resultpath = "result//212007_flat+9.fit"
fits.writeto(resultpath, data_stack, header=header1, overwrite=True)


