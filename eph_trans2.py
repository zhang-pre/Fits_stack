from scipy import interpolate as inter

f_header = open("0209Header.txt", "r+")
f_eph = open("0209az.txt", "r+")
f_result = open("result_20210209_az.txt", "w+")

f_headerlist = f_header.readlines()
f_ephlist = f_eph.readlines()

for header in f_headerlist:
    header = header.split()
    headertime = header[2].split(":")
    filename = header[0]
    headertotal = int(headertime[0]) * 3600 + int(headertime[1]) * 60 + int(headertime[2]) + 10.0
    print(headertotal)

    pre1_ra = 0
    pre1_de = 0
    pre2_ra = 0
    pre2_de = 0
    post1_ra = 0
    post1_de = 0
    post2_ra = 0
    post2_de = 0
    pre1time = 0
    pre2time = 0
    post1time = 0
    post2time = 0

    for ephline in f_ephlist:
        ephline = ephline.split()
        ephtimetotal = int(ephline[3]) * 3600 + int(ephline[4]) * 60
        # ra_h = int(ephline[5])
        # ra_m = int(ephline[6])/60.0
        # ra_s = float(ephline[7])/3600.0
        #
        # de_d = int(ephline[8])
        # de_m = int(ephline[9])/60.0
        # de_s = float(ephline[10])/3600.0

        ra = float(ephline[5])
        de = float(ephline[6])

        hme = headertotal - ephtimetotal

        if hme >= 0 and hme < 60:
            pre1_ra = ra
            pre1_de = de
            pre1time = ephtimetotal
        if hme >= 60 and hme <= 120:
            pre2_ra = ra
            pre2_de = de
            pre2time = ephtimetotal
        if hme >= (-60.0) and hme < 0:
            post1_ra = ra
            post1_de = de
            post1time = ephtimetotal
        if hme >= (-120.0) and hme < (-60):
            post2_ra = ra
            post2_de = de
            post2time = ephtimetotal

    ralist = [pre2_ra, pre1_ra, post1_ra, post2_ra]
    delist = [pre2_de, pre1_de, post1_de, post2_de]
    timelist = [pre2time, pre1time, post1time, post2time]
    funra = inter.interp1d(timelist, ralist, kind="quadratic")
    ra_inter = funra(headertotal)
    funde = inter.interp1d(timelist, delist, kind="quadratic")
    de_inter = funde(headertotal)

    f_result.write(filename + "   " + str(ra_inter) + "   " + str(de_inter) + "   " + "\n")