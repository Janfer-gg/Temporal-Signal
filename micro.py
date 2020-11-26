
import rpy2.robjects as robjects
import MySQLdb

gene='thrl'
species='Escherichia coli BL21(DE3)'
tags = ''
ID=''
db = MySQLdb.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8' )
cursor = db.cursor()
# sql = "SELECT * FROM micro WHERE gene = %s AND species = 'K-12 substr. MG1655' " ,gene

# 执行SQL语句
cursor.execute("SELECT * FROM micro WHERE species = '%s'  AND (gene = '%s' or locus_tag = '%s' or ID = '%s')  " %(species,gene,tags,ID))

# 获取所有记录列表
results = cursor.fetchone()

print(results)
gene = results[1]
# length = int(results[3]) - int(results[2])
# seq = results[5]
# GC = results[6]
# h_start = results[7].split(";")
# h_end = results[8].split(";")
gRNA1 = results[10]
# strand1 = results[10]
# pos1 = results[11]
gRNA2 = results[14]
# strand2 = results[13]
# pos2 = results[14]
# destroy = results[16]
# up = results[17]
# up_dir = results[18]
# up_dis = results[19]
# down = results[20]
# down_dir = results[21]
# down_dis = results[22]
# 打印结果
print(gRNA1, gRNA2)


# 关闭数据库连接
db.close()

# filepath = 'C://Users//41518//Desktop'
#
# robjects.r(r"source('C://Users//41518//Desktop//work//Ubigene//micro-organism//micro_plot.R')")
# robjects.r['micro_plot'](gene,length,seq,filepath,pos1,pos2)



# reader = csv.reader(open('C://Users//41518//Desktop//微生物//micro.csv','r'))
# list_result = list(reader)
#
# species='K-12 substr. MG1655'
# gene='carB'
# filepath = 'C://Users//41518//Desktop//微生物//'+gene
#
#
# for results in list_result:
#     if results[0]==species and results[1]==gene :
#         length = int(results[3])-int(results[2])
#         GC = results[6]
#         h_start = results[7]
#         h_end = results[8]
#         gRNA1 = results[9]
#         strand1 = results[10]
#         pos1 = results[11]
#         gRNA2 = results[12]
#         strand2 = results[13]
#         pos2 = results[14]
#
# print(length,GC,gRNA1,strand1,pos1,gRNA2,strand2,pos2,h_start,h_end)
#

