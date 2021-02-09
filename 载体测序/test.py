import os
import csv
import json

# os.system(r'D:\R\R-4.0.1\bin\x64\Rscript ./test.R')

with open('result.csv', 'r', encoding='utf8')as aa:
    reader = csv.reader(aa)
    table = []
    for line in reader:
        table.append(line)
    reuslt=json.dumps(table[1:])
    print(reuslt)

