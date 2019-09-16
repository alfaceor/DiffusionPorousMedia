#!/usr/bin/python
# Read files from csv file and check if file exists?
import csv
from SimulData import *

with open("RunJobs.csv") as csvfile:
  mpg = list(csv.DictReader(csvfile))

for i in range(len(mpg)):
  sd = SimulData()
  sd.metadataFromCSVentry("", mpg[i])
  sd.defineFiles()
  print sd.arguments


