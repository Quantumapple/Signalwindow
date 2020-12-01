#!/usr/bin/python

import os

InputDir = "/xrootd/store/user/jongho/Delphes/Bkg/0009"
os.system("ls " + InputDir + " > ./inputlist.txt")

myFile = open('./inputlist.txt','r')
lines = myFile.readlines()
newFile = open('./inputlist_my.txt','w')
for name in lines:
    print name
    newFile.write("root://cms-xrdr.private.lo:2094///xrd/store/user/jongho/Delphes/Bkg/0009/" + name)

myFile.close()
newFile.close()


#line_count = len(myFile.readlines())
