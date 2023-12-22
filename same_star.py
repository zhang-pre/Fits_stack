
f1 = open("20210415_AM.txt")
f2 = open("20210415_our.txt")

f3 = open("20210415_AM_clip.txt", mode="w+")
f4 = open("20210415_our_clip.txt", mode="w+")

f1list = f1.readlines()
f2list = f2.readlines()

for line1 in f1list:

    starno1 = line1.split()[0]

    for line2 in f2list:

        starno2 = line2.split()[0]

        if starno1 == starno2:

            f3.write(line1)
            f4.write(line2)

f1.close()
f2.close()
f3.close()
f4.close()