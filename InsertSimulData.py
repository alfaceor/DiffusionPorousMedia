#!/usr/bin/python

import sqlite3

db = sqlite3.connect("SimulData.db")

cursor = db.cursor()

Dir = "RandomIniCond_02/"

strF2BList  = [ "1.0", "0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001" ] 
strB2FList  = [ "1.0", "0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001" ]
strB2BList  = [ "1.0", "0.1", "0.01", "0.001", "0.0001", "0.00001", "0.000001" ]
HradList = [ 10, 20, 40, 80, 160, 320  ]
BradList = [ i+1 for i in HradList ]
strHradList = map(str, HradList)
strBradList = map(str, BradList)

count = 0
for hr in range(len(strHradList)):
  strHrad = strHradList[hr]
  strBrad = strBradList[hr]
  for f2b in range(len(strF2BList)):
    strF2B  = strF2BList[f2b]
    for b2f in range(len(strB2FList)):
      strB2F  = strB2FList[b2f]
      for b2b in range(len(strB2BList)):
        strB2B  = strB2BList[b2b]
        print Dir, strF2B, strB2F, strB2B, strHrad, strBrad, count
        count = count +1
        cursor.execute('''INSERT INTO Simulation(Directory, strF2B, strB2F, strB2B, strHrad, strBrad) VALUES(?, ?, ?, ?, ?, ?)''', (Dir, strF2B, strB2F, strB2B, strHrad, strBrad))

db.commit()
db.close()
