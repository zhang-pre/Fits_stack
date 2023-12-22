
from constants_calculate import *
from astropy.io import fits
from test import *



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

x_list_1 = np.array([939.059, 973.342, 571.479])
y_list_1 = np.array([176.852, 278.677, 520.456])

x_list_2 = np.array([939.968, 974.246, 572.955])
y_list_2 = np.array([193.104, 294.906, 535.507])

x_list_3 = np.array([931.979, 966.276, 560.433])
y_list_3 = np.array([209.06, 310.832, 551.761])

x_list_4 = np.array([935.418, 969.76, 568.349])
y_list_4 = np.array([209.06, 326.012, 566.751])

x_list_5 = np.array([934.097, 968.509, 567.108])
y_list_5 = np.array([241.021, 342.739, 583.493])

x_list_6 = np.array([914.057, 948.363, 547.088])
y_list_6 = np.array([257.435, 359.114, 599.917])

x_list_7 = np.array([904.702, 939.023, 537.698])
y_list_7 = np.array([273.779, 375.506, 616.279])

# 其他图像转换到第一幅图像
p0 = [0, 0, 0]
Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_2))
a, b, c = Para_ksi[0]
Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_2))
e, f, g = Para_ita[0]
const_six = [a, b, c, e, f, g]

print(const_six)

p0 = [0, 0, 0]
Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_3))
a, b, c = Para_ksi[0]
Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_3))
e, f, g = Para_ita[0]
const_six = [a, b, c, e, f, g]

print(const_six)

p0 = [0, 0, 0]
Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_4))
a, b, c = Para_ksi[0]
Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_4))
e, f, g = Para_ita[0]
const_six = [a, b, c, e, f, g]

print(const_six)

p0 = [0, 0, 0]
Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_5))
a, b, c = Para_ksi[0]
Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_5))
e, f, g = Para_ita[0]
const_six = [a, b, c, e, f, g]

print(const_six)

p0 = [0, 0, 0]
Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_6))
a, b, c = Para_ksi[0]
Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_6))
e, f, g = Para_ita[0]
const_six = [a, b, c, e, f, g]

print(const_six)

p0 = [0, 0, 0]
Para_ksi = leastsq(error_ksi, p0, args=(x_list_1, y_list_1, x_list_7))
a, b, c = Para_ksi[0]
Para_ita = leastsq(error_ita, p0, args=(x_list_1, y_list_1, y_list_7))
e, f, g = Para_ita[0]
const_six = [a, b, c, e, f, g]

print(const_six)