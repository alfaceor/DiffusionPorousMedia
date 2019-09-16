#!/usr/bin/python

import sqlite3

db = sqlite3.connect("SimulData.db")

cursor = db.cursor()

Dir = "RandomIniCond_02/"

strF2B  = "1.0"
strB2F  = "1.0"
strB2B  = "1.0"
strHrad = "10"
strBrad = "11"


aaa = cursor.execute('''SELECT FitMin, FitMax FROM Simulation WHERE Directory = ? AND strF2B = ? AND strB2F = ? AND strB2B = ? AND strHrad = ? AND strBrad =?''', (Dir, strF2B, strB2F, strB2B, strHrad, strBrad))

print aaa
#print aaa.fetchall()
FitMin, FitMax = aaa.fetchone()
print FitMin, FitMax
#for row in aaa.fetchall():
#  print row, "aaa"

#db.commit()
db.close()



