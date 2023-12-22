import os

import numpy as np
from astropy.coordinates import Angle
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.io import fits

# fwhm = 5.5

#todo 读取positions目录下的star-position.txt文件，把里面的前三列转换成赤经，后三列转换成赤纬，并存储
with open('positions/5.27_S9.txt', 'r') as file:
    lines = file.readlines()

# 解析每行数据并将其转换为赤经和赤纬
ra_deg = []
dec_deg = []
for line in lines:
    values = line.strip().split(',')
    ra_str , dec_str= values[0],values[1]
    # 使用Angle对象将字符串转换为度
    ra_deg.append(Angle(ra_str).degree)
    dec_deg.append(Angle(dec_str).degree)

# 现在ra_deg和dec_deg分别包含了赤经和赤纬的度数值
#todo 读取wcs目录下的fits文件，输入前面存储的赤经赤纬，wcs转换成像素坐标
my_filename = "wcs/5.27/"
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
filelist = walkFile(my_filename)



# todo  偏移的话，应该用background填补;如果让偏移量是浮点数，并且运用插值的方法，是不是更准确点呢？
# 读取第一个FITS文件来获取WCS信息
first_image = fits.open(filelist[0])
first_wcs = WCS(first_image[0].header)
ori_x,ori_y = first_wcs.all_world2pix(ra_deg[0], dec_deg[0], 0)
# 初始化一个累积图像
accumulated_image = np.zeros(first_image[0].data.shape)

# 遍历所有FITS文件
for i, filename in enumerate(filelist):
    hdul = fits.open(filename)
    cur_data = hdul[0].data
    wcs = WCS(hdul[0].header)
    mean, median, std = sigma_clipped_stats(cur_data, sigma=3.0)
    background = 3 * median - 2 * mean
    # 计算像素坐标的偏移
    if i > 0:
        x,y = wcs.all_world2pix(ra_deg[i], dec_deg[i], 0)
        dx, dy = int(ori_x-x), int(ori_y-y)
    else:
        dx,dy = 0,0

    # 加上像素坐标偏移量并叠加图像,
    # 先加上偏移量，再差值拟合，再叠加
    # offset_data = np.roll(cur_data, (dx,dy), axis=(0, 1))
    offset_data = np.roll(cur_data, (dx, dy), axis=(0, 1))
    # 在进行偏移操作后，让越界的数据填充为0而不是从尾又移到头
    if dx > 0:
        offset_data[0:dx, :] = background
    elif dx < 0:
        offset_data[dx:, :] = background
    if dy > 0:
        offset_data[:, :dy] = background
    elif dy < 0:
        offset_data[:, dy:] = background

    accumulated_image += offset_data

    hdul.close()
# 保存叠加后的图像

output_filename = "accumulated_image27_1+5_2.fits"
tempfilename = "original/20230526S9-I_001.fits"
newHeader = fits.open(tempfilename)[0].header

fits.writeto(output_filename, accumulated_image, header=newHeader, overwrite=True)

# 关闭第一个FITS文件
first_image.close()
#todo 后面的每一个像素坐标和第一个像素坐标相减，得到像素坐标的偏移量
#todo 后面的fits图片累加到第一张图片上，但是要注意加上像素坐标偏移量
#todo 对坐标插值拟合

