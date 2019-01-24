"""
Code to look at Voyager 2 Plasma Wave Subsystem data around the Neptune ring plane crossing,
to compare it with that broadcast on Neptune All Night
"""

import csv
import datetime
import numpy as np
import matplotlib.pyplot as plt

img=[]
tstart=datetime.datetime.strptime("1989-237T02:00:00.000000","%Y-%jT%H:%M:%S.%f")
tstop =datetime.datetime.strptime("1989-237T03:00:00.000000","%Y-%jT%H:%M:%S.%f")
with open('../../Data/Voyager/pws-neptune/PWS_NEP_SA_4SEC.TAB', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|',skipinitialspace=True)
    for row in spamreader:
        time=datetime.datetime.strptime(row[0]+"000","%Y-%jT%H:%M:%S.%f")
        if time>=tstart and time<tstop:
            this_row=[]
            for channel in row[2:]:
                this_row.append(int(channel))
            print(time)
            img.append(this_row)
img=np.array(img)
plt.imshow(img.transpose())
plt.show()

