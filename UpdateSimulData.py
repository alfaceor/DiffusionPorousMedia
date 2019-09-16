#!/usr/bin/python

import sqlite3

db = sqlite3.connect("SimulData.db")

cursor = db.cursor()

Dir = "RandomIniCond_02/"

HradList = [ 10, 20, 40, 80, 160, 320 ] 
BradList = [ i+1 for i in HradList ] 
strHradList = map(str, HradList) 
strBradList = map(str, BradList) 
strF2BList  = [ "1.0" ]
strB2FList  = [ "1.0", "0.1", "0.01", "0.001", "0.001", "0.0001", "0.00001", "0.000001"   ] 
strB2BList  = [ "1.0", "0.1", "0.01", "0.001", "0.001", "0.0001", "0.00001", "0.000001"   ] 

FitMin  = 100000
FitMax  = 1000000

for hr  in range(len(strHradList)):
  for f2b in range(len(strF2BList)):
    for b2f in range(len(strB2FList)):
      for b2b in range(len(strB2BList)):
        strF2B  = strF2BList[f2b]  
        strB2F  = strB2FList[b2f]  
        strB2B  = strB2BList[b2b]  
        strHrad = strHradList[hr] 
        strBrad = strBradList[hr] 
        cursor.execute('''UPDATE Simulation SET FitMin = ?, FitMax = ? WHERE Directory = ? AND strF2B = ? AND strB2F = ? AND strB2B = ? AND strHrad = ? AND strBrad =?''', (FitMin, FitMax, Dir, strF2B, strB2F, strB2B, strHrad, strBrad))

db.commit()
db.close()



