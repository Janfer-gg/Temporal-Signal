#!/usr/bin/env python
#-*- coding:utf-8 -*-
# author:NAIVE
# datetime:2020/12/25 14:41
import weasyprint
from html.parser import HTMLParser
import os
import django
import os, sys

sys.path.insert(0, r'F:\pythoncode\auto_plan_system\AutoPlanSys')
print(sys.path)
import AutoPlanSys
import pdfkit
import xlrd

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "AutoPlanSys.settings")
django.setup()
# config = pdfkit.configuration(wkhtmltopdf=r"E:\tools\wkhtmltopdf\bin\wkhtmltopdf.exe")
from django.template.loader import render_to_string

import pymysql

conn = pymysql.connect(host="47.92.215.36",
                            user="root",
                            password="@ubigene&yuanjing@",
                            db="ubigene_test",
                            port=3306,
                            charset="utf8")
cursor = conn.cursor()
def admin_order_pdf(file_path, content):
    html_parser = HTMLParser()
    new_html = html_parser.unescape(content)
    weasyprint.HTML(string=new_html, base_url=r'E:\workspace\add_pic\gene_pdf_pac').write_pdf(file_path)


def gene_pdf_item(file_path,plasmid_name,gene_name,seqs,item_num,gene_id,chunk_all_grna_list,ri_list):

    seq_list = [seq.replace('(','[').replace(')',']') for seq in seqs]


    item_seq_list = zip(item_num,seq_list)
    zip_chunk_all_grna_list = zip(chunk_all_grna_list,ri_list)

    html = render_to_string(r'synthego/nav_contract.html', {"plasmid_name": plasmid_name,"item_seq_list": item_seq_list,"gene_name": gene_name,'gene_id':gene_id,'zip_chunk_all_grna_list':zip_chunk_all_grna_list,'ri_list':ri_list})
    admin_order_pdf(file_path, html)

    with open('test_html.html','w',encoding='utf-8')as fff:
        fff.write(html)


with open('../gene_name_id')as ff:
    ff_list = ff.readlines()

with open('YKO.txt','r')as fY:
    YKO_str = fY.read().strip().upper()


none_c1 = YKO_str[0:14]
U6_promoter = YKO_str[14:263]
gRNA = YKO_str[263:284]
scaffold = YKO_str[284:360]
none_c2 = YKO_str[360:372]
cbh_promoter = YKO_str[372:1170]
none_c3 = YKO_str[1170:1182]
spcas9 = YKO_str[1182:5454]
none_c4 = YKO_str[5454:5484]
bGH_poly = YKO_str[5484:5692]
none_c5 = YKO_str[5692:5721]
cmv_promoter = YKO_str[5721:6310]
none_c6 = YKO_str[6310:6348]
EGFP_PA_PR = YKO_str[6348:7731]
none_c7 = YKO_str[7731:7746]
sv40_ply_singal = YKO_str[7746:7868]
none_c8 = YKO_str[7868:9022]
AMPR = YKO_str[9022:9883]
none_c9 = YKO_str[9883:10053]
ORI = YKO_str[10053:10642]


all_grna_list = []

for none_1 in none_c1:
    grna_color_dict = {}
    grna_color_dict[none_1]= 'none'
    all_grna_list.append(grna_color_dict)
for U6_p in U6_promoter:
    grna_U6_p_dict = {}
    grna_U6_p_dict[U6_p]= '#523cf0'
    all_grna_list.append(grna_U6_p_dict)

print('-'*100)
print(gRNA)
for gc, gR in enumerate(gRNA):
    grna_gR_dict = {}
    if gc == 0 :
        grna_gR_dict['g'] = '#b2bac7'
    else:
        grna_gR_dict['N']= '#b2bac7'
    all_grna_list.append(grna_gR_dict)



for scaf in scaffold:
    grna_scaf_dict = {}
    grna_scaf_dict[scaf] = '#b2bac7'
    all_grna_list.append(grna_scaf_dict)

for none_2 in none_c2:
    grna_none_2_dict = {}
    grna_none_2_dict[none_2] = 'none'
    all_grna_list.append(grna_none_2_dict)

for cbh_p in cbh_promoter:
    grna_cbh_p_dict = {}
    grna_cbh_p_dict[cbh_p] = '#c18ece'
    all_grna_list.append(grna_cbh_p_dict)
for none_3 in none_c3:
    grna_none_3_dict = {}
    grna_none_3_dict[none_3] = 'none'
    all_grna_list.append(grna_none_3_dict)
for spc in spcas9:
    grna_spc_dict = {}
    grna_spc_dict[spc] = '#ff6325'
    all_grna_list.append(grna_spc_dict)
for none_4 in none_c4:
    grna_none_4_dict = {}
    grna_none_4_dict[none_4] = 'none'
    all_grna_list.append(grna_none_4_dict)
for bGH_p in bGH_poly:
    grna_bGH_p_dict = {}
    grna_bGH_p_dict[bGH_p] = '49eff5'
    all_grna_list.append(grna_bGH_p_dict)
for none_5 in none_c5:
    grna_none_5_dict = {}
    grna_none_5_dict[none_5] = 'none'
    all_grna_list.append(grna_none_5_dict)
for cmv_p in cmv_promoter:
    grna_cmv_p_dict = {}
    grna_cmv_p_dict[cmv_p] = '#044a7d'
    all_grna_list.append(grna_cmv_p_dict)
for none_6 in none_c6:
    grna_none_6_dict = {}
    grna_none_6_dict[none_6] = 'none'
    all_grna_list.append(grna_none_6_dict)
for EGFP_P in EGFP_PA_PR:
    grna_EGFP_P_dict = {}
    grna_EGFP_P_dict[EGFP_P] = '#9c7d69'
    all_grna_list.append(grna_EGFP_P_dict)
for none_7 in none_c7:
    grna_none_7_dict = {}
    grna_none_7_dict[none_7] = 'none'
    all_grna_list.append(grna_none_7_dict)
for sv40_p in sv40_ply_singal:
    grna_sv40_p_dict = {}
    grna_sv40_p_dict[sv40_p] = '#30a47d'
    all_grna_list.append(grna_sv40_p_dict)
for none_8 in none_c8:
    grna_none_8_dict = {}
    grna_none_8_dict[none_8] = 'none'
    all_grna_list.append(grna_none_8_dict)
for AM_P in AMPR:
    grna_AM_P_dict = {}
    grna_AM_P_dict[AM_P] = '#eaa62f'
    all_grna_list.append(grna_AM_P_dict)
for none_9 in none_c9:
    grna_none_9_dict = {}
    grna_none_9_dict[none_9] = 'none'
    all_grna_list.append(grna_none_9_dict)
for ORI_p in ORI:
    grna_ORI_p_dict = {}
    grna_ORI_p_dict[ORI_p] = '#ffe056'
    all_grna_list.append(grna_ORI_p_dict)

print(all_grna_list)


def chunkIt(seq, avg):
    avg = avg
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


print(chunkIt(all_grna_list, 10))
chunk_all_grna_list = chunkIt(chunkIt(all_grna_list, 10), 10)
print(len(chunk_all_grna_list))

ri_list = [ri for ri in range(1,10700,100)]
print(len(chunk_all_grna_list))
import time

for ff in ff_list:
    s = time.time()

    ff_str = ff.strip()
    gene_name_id_list = ff_str.split('	')

    print(gene_name_id_list)

    plasmid_name = gene_name_id_list[0]
    gene_id = gene_name_id_list[1]
    gene_name = gene_name_id_list[2]

    cursor.execute('select id from dede_archives where title like "%{}%"'.format(plasmid_name))
    id_tuple = cursor.fetchone()
    idd = id_tuple[0]
    cursor.execute('select itemno,serialnumber from dede_addon21  where aid="{}"'.format(idd))
    itemno, serialnumber = cursor.fetchone()
    print(gene_name,itemno,serialnumber)

    # serialnumber_list = serialnumber.split()
    # ser_seq_20_list = []
    # for serialnum in serialnumber_list:
    #     ser_seq = serialnum.split('(')[0]
    #
    #     ser_seq_20 = ser_seq[0:20]
    #     ser_seq_20_list.append(ser_seq_20)



    file_path = rf'E:\workspace\add_pic\gene_pdf_pac\detail_pdf_gRNA\{plasmid_name}.pdf'

    gene_pdf_item(file_path,plasmid_name,gene_name,serialnumber.split(),itemno.split(),gene_id,chunk_all_grna_list,ri_list)

    end = time.time()


    t = end-s

    print(t)