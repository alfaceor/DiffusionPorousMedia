#!/usr/bin/python

strProgram="./rw3D_RegularPorous__time_1e7__RandomIniCond_02 "
strF2B = "1.0"
strB2F = "0.001"
strB2B = "0.000001"
HradList = [20, 40 ] #, 80, 160, 320]
for hr in range(len(HradList)):
  strHrad = str(HradList[hr])
  strBrad = str(HradList[hr]+1)
  print strProgram+" --pF2B "+strF2B+" --pB2F "+strB2F+" --pB2B "+strB2B+" --Hradius "+strHrad+ " --Bradius "+strBrad
