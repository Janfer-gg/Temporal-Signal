# import rpy2.robjects as robjects
# import csv
# robjects.r(r"source('C://Users//41518//Desktop//work//Ubigene//primer//primer_tool.R')")

# filepath ='C://Users//41518//Desktop//靶位点4/SAMD8'
# term = 'SAMD8'
# species = 'Human'
# gRNA1 = 'AGATGGGAGTCTACCTGAAG'
# gRNA2 = ''
#
# robjects.r['primer_tool'](term,species,gRNA1,gRNA2,filepath)
#
# result = csv.reader(open('{}/result.csv'.format(filepath),'r'))
# list_result = list(result)
#
# result_dict = { i:j for i,j in zip(list_result[0],list_result[1])}
# print(result_dict)

#!/usr/bin/env python
#-*- coding:utf-8 -*-
# author:NAIVE
# datetime:2020/10/30 12:00
import os
import csv
#
filepath ='C://Users//41518//Desktop//靶位点4/SAMD8'
term = 'SAMD8'
species = 'Human'
gRNA1 = 'AGATGGGAGTCTACCTGAAG'
gRNA2 = 'NULL'
work_path = 'C://Users//41518//Desktop//work//Ubigene//primer'
os.system('D://R//R-4.0.1//bin//x64//Rscript C://Users//41518//Desktop//primer_tool.R'+ " " + term + " " + species + " " + gRNA1 + " " + gRNA2 + " " + filepath+ " " + work_path)

result = csv.reader(open('{}/result.csv'.format(filepath),'r'))
list_result = list(result)

result_dict = { i:j for i,j in zip(list_result[0],list_result[1])}
print(result_dict)