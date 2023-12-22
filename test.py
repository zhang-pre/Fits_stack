import numpy as np
import os
# 双线性插值
def merge_plate_6model(fit_data, consts, x_length, y_length):

    data = np.zeros(x_length * y_length).reshape(x_length, y_length)
    # b = [([0] * len(fit_data[0])) for i in range(len(fit_data))]
    yy = len(fit_data)
    xx = len(fit_data[0])
    i = 0
    while i < xx:
        j = 0
        while j < yy:
            x = consts[0] * i + consts[1] * j + consts[2]
            y = consts[3] * i + consts[4] * j + consts[5]
            if x >= 0 and y >= 0:
                lx = int(x)
                ly = int(y)
                if lx >= 0 and ly >= 0 and lx + 1 < xx and ly + 1 < yy:
                    xd = x - lx
                    yd = y - ly
                    fe = 0
                    ff = 0
                    a1 = 1 - xd
                    a2 = lx + 1
                    a3 = ly + 1
                    fe = a1 * fit_data[ly][lx] + xd * fit_data[ly][a2]
                    ff = a1 * fit_data[a3][lx] + xd * fit_data[a3][a2]
                    # b[j][i] = (1 - yd) * fe + yd * ff
                    data[j][i] = (1 - yd) * fe + yd * ff
            j = j + 1
        i = i + 1
    return data

def merge_plate_12model(fit_data, consts, x_length, y_length):

    data = np.zeros(x_length * y_length).reshape(x_length, y_length)
    # b = [([0] * len(fit_data[0])) for i in range(len(fit_data))]
    yy = len(fit_data)
    xx = len(fit_data[0])
    i = 0
    while i < xx:
        j = 0
        while j < yy:
            x = consts[0] * i * i + consts[1] * j * j + consts[2] * i * j + consts[3] * i + consts[4] * j + consts[5]
            y = consts[6] * i * i + consts[7] * j * j + consts[8] * i * j + consts[9] * i + consts[10] * j + consts[11]
            if x >= 0 and y >= 0:
                lx = int(x)
                ly = int(y)
                if lx >= 0 and ly >= 0 and lx + 1 < xx and ly + 1 < yy:
                    xd = x - lx
                    yd = y - ly
                    fe = 0
                    ff = 0
                    a1 = 1 - xd
                    a2 = lx + 1
                    a3 = ly + 1
                    fe = a1 * fit_data[ly][lx] + xd * fit_data[ly][a2]
                    ff = a1 * fit_data[a3][lx] + xd * fit_data[a3][a2]
                    # b[j][i] = (1 - yd) * fe + yd * ff
                    data[j][i] = (1 - yd) * fe + yd * ff
            j = j + 1
        i = i + 1
    return data

def merge_plate_20model(fit_data, consts, x_length, y_length):

    data = np.zeros(x_length * y_length).reshape(x_length, y_length)
    # b = [([0] * len(fit_data[0])) for i in range(len(fit_data))]
    yy = len(fit_data)
    xx = len(fit_data[0])
    i = 0
    while i < xx:
        j = 0
        while j < yy:
            x = consts[0] * i * i * i + consts[1] * j * j * j + consts[2] * i * i * j + consts[3] * i * j * j + \
                consts[4] * i * i + consts[5] * j * j + consts[6] * i * j + consts[7] * i + consts[8] * j + consts[9]

            y = consts[10] * i * i * i + consts[11] * j * j * j + consts[12] * i * i * j + consts[13] * i * j * j + \
                consts[14] * i * i + consts[15] * j * j + consts[16] * i * j + consts[17] * i + consts[18] * j + consts[19]

            if x >= 0 and y >= 0:
                lx = int(x)
                ly = int(y)
                if lx >= 0 and ly >= 0 and lx + 1 < xx and ly + 1 < yy:
                    xd = x - lx
                    yd = y - ly
                    fe = 0
                    ff = 0
                    a1 = 1 - xd
                    a2 = lx + 1
                    a3 = ly + 1
                    fe = a1 * fit_data[ly][lx] + xd * fit_data[ly][a2]
                    ff = a1 * fit_data[a3][lx] + xd * fit_data[a3][a2]
                    # b[j][i] = (1 - yd) * fe + yd * ff
                    data[j][i] = (1 - yd) * fe + yd * ff
            j = j + 1
        i = i + 1
    return data

def merge_plate_30model(fit_data, consts, x_length, y_length):

    data = np.zeros(x_length * y_length).reshape(x_length, y_length)
    # b = [([0] * len(fit_data[0])) for i in range(len(fit_data))]
    yy = len(fit_data)
    xx = len(fit_data[0])
    i = 0
    while i < xx:
        j = 0
        while j < yy:
            x = consts[0] * i * i * i * i + consts[1] * j * j * j * j + consts[2] \
                * i * i * i * j + consts[3] * i * j * j * j + consts[4] * i * i * j * j +\
                consts[5] * i * i * i + consts[6] * j * j * j + consts[7] * i * i * j + \
                consts[8] * i * j * j + consts[9] * i * i + consts[10] * j * j + consts[11] * i * j + \
                consts[12] * i + consts[13] * j + consts[14]

            y = consts[15] * i * i * i * i + consts[16] * j * j * j * j + consts[17] \
                * i * i * i * j + consts[18] * i * j * j * j + consts[19] * i * i * j * j +\
                consts[20] * i * i * i + consts[21] * j * j * j + consts[22] * i * i * j + \
                consts[23] * i * j * j + consts[24] * i * i + consts[25] * j * j + consts[26] * i * j + \
                consts[27] * i + consts[28] * j + consts[29]

            if x >= 0 and y >= 0:
                lx = int(x)
                ly = int(y)
                if lx >= 0 and ly >= 0 and lx + 1 < xx and ly + 1 < yy:
                    xd = x - lx
                    yd = y - ly
                    fe = 0
                    ff = 0
                    a1 = 1 - xd
                    a2 = lx + 1
                    a3 = ly + 1
                    fe = a1 * fit_data[ly][lx] + xd * fit_data[ly][a2]
                    ff = a1 * fit_data[a3][lx] + xd * fit_data[a3][a2]
                    # b[j][i] = (1 - yd) * fe + yd * ff
                    data[j][i] = (1 - yd) * fe + yd * ff
            j = j + 1
        i = i + 1
    return data



filename = "original/"

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