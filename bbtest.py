import numpy as np

# 假设你有一个二维数组 data，和 x、y 方向的偏移量 dx 和 dy
data = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9]])
dx = -2  # x 方向的偏移量
dy = -1  # y 方向的偏移量

# 使用 np.roll 函数来进行偏移

offset_data = np.roll(data, (dx, dy), axis=(0, 1))
# 在进行偏移操作后，让越界的数据填充为0而不是从尾又移到头
if dx>0:
    offset_data[0:dx,:] = 0
elif dx<0:
    offset_data[dx:,:] = 0
if dy > 0:
    offset_data[:, :dy] = 0
elif dy < 0:
    offset_data[:, dy:] = 0
print("原始数组:")
print(data)
print("偏移后的数组:")
print(offset_data)



# import numpy as np
#
# # 创建一个示例数组
# arr = np.array([[1, 2, 3],
#                 [4, 5, 6],
#                 [7, 8, 9]])
#
# # 定义偏移量
# x,y =1,1
#
# # 计算要填充的0的数量
# fill_count = abs(shift) % len(arr)
#
# # 计算正向或负向偏移
# if shift < 0:
#     # 负向偏移，数据向左滚动
#     result = np.roll(arr, shift)
#     result[fill_count:] = 0  # 填充后面的位置为0
# else:
#     # 正向偏移，数据向右滚动
#     result = np.roll(arr, shift)
#     result[:fill_count] = 0  # 填充前面的位置为0
#
# print(result)
import os

import numpy as np
from astropy.coordinates import Angle
from astropy.stats import sigma_clipped_stats
# from astropy.wcs import WCS
# from astropy.io import fits
#
#  # 打开 FITS 文件
# hdulist = fits.open('result_star/20230526S9-I_001_6+4.fit')
# data = hdulist[0].data
# print(data.shape)
# # 获取 WCS 信息
# wcs = WCS(hdulist[0].header)  # 这里假设 WCS 信息存储在 FITS 文件的头文件中
#
# # 定义天球坐标（赤道坐标），例如一个天体的赤经和赤纬
# # 339.27902211, -10.39047485
# ra = 339.2791619999999  # 以度为单位的赤经
# dec = -10.390609723333332  # 以度为单位的赤纬
# ra2 = 339.2792051916666
# dec2 = -10.390598058611111
# # 将天球坐标转换为像素坐标
# x, y = wcs.all_world2pix(ra, dec, 0)  # 第三个参数通常是 0，表示使用最精确的坐标映射
# x2, y2 = wcs.all_world2pix(ra2, dec2, 0)  # 第三个参数通常是 0，表示使用最精确的坐标映射
# # x 和 y 现在包含了对应的像素坐标
# print(f'赤经 {ra}, 赤纬 {dec} 对应的像素坐标为 x={x}, y={y}')
# print(f'赤经 {ra2}, 赤纬 {dec2} 对应的像素坐标为 x={x2}, y={y2}')

