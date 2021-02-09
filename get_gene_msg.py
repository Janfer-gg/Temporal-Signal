#!/usr/bin/env python
#-*- coding:utf-8 -*-
# author:NAIVE
# datetime:2020/9/17 16:59
import requests
from scrapy.selector import Selector
import os

import logging

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)s %(levelname)s %(message)s",
                    datefmt='%Y-%m-%d  %H:%M:%S %a')

term_get_id_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}&sort=Relevance&retmode=xml'

id_get_esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&retmode=xml'

img_url = 'C://Users//41518//Desktop//work//Ubigene//gene_img//{}'
plan_url = "C://Users//41518//Desktop//work//Ubigene//R//"
run_path = "C://Users//41518//Desktop//work//Ubigene//R"
file_url = 'C://Users//41518//Desktop//work//Ubigene//gene_img//{}//'
from requests.adapters import HTTPAdapter
s = requests.Session()
requests.packages.urllib3.disable_warnings()
s.mount('http://', HTTPAdapter(max_retries=100))
s.mount('https://', HTTPAdapter(max_retries=100))

def get_resp_text(url, **kwargs):
    resp = s.get(url=url,verify=False)
    return resp.text


def xml_selector(api_url):
    xml = get_resp_text(api_url)
    sel = Selector(type='xml', text=xml)

    return sel


def get_id(term):
    sel = xml_selector(api_url=term_get_id_url.format(term))
    id_list = sel.xpath('//Id//text()').extract()

    return id_list

with open(r'C:\Users\41518\Desktop\work\Ubigene\abcam','r')as f:

    abcam_gene_list = f.readlines()

for abcam_g in abcam_gene_list:

    abcam_gene = abcam_g.strip()

    try:

        idd = get_id(abcam_g)[0]
    except:
        continue
    # sel_esearch = xml_selector(api_url=id_get_esearch_url.format(idd))
    #
    # # esearch_text = get_resp_text(url=id_get_esearch_url.format(idd))
    # # esearch_res = re.search(r'type (.*?),',esearch_text).group()
    # # item['GeneType'] = esearch_res.replace('type ','').replace(',','')
    # EnsemblIdList = sel_esearch.xpath('//Entrezgene_gene//Object-id_str//text()').extract()
    # for i in EnsemblIdList:
    #     if i.startswith('EN') and len(i) > 9:
    #         EnsemblId = i
    #         break
    #     else:
    #         EnsemblId = ''
    #
    # print(idd,EnsemblId)
    # planA = '{}planA.R'.format(plan_url)
    # planB = '{}PlanB.R'.format(plan_url)
    # planC = '{}PlanC.R'.format(plan_url)
    # image_file_path = img_url.format(idd)
    # isexists = os.path.exists(image_file_path)
    # if not isexists:
    #     os.mkdir(image_file_path)
    # else:
    #     continue
    # os.system(r'Rscript {}'.format(planA) + " " + 'Human' + " " + EnsemblId + " " + image_file_path + " " + run_path + " " + abcam_g)
    # os.system(r'Rscript {}'.format(planB) + " " + 'Human' + " " + EnsemblId + " " + image_file_path + " " + run_path + " " + abcam_g)
    # os.system(r'Rscript {}'.format(planC) + " " + 'Human' + " " + EnsemblId + " " + image_file_path + " " + run_path + " " + abcam_g)
    #
    # is_result_path_a10 = '{}A10%.txt'.format(file_url.format(idd))
    # is_result_path_a20 = '{}A20%.txt'.format(file_url.format(idd))
    # is_result_path_a60 = '{}A60%.txt'.format(file_url.format(idd))
    # is_result_path_a90 = '{}A90%.txt'.format(file_url.format(idd))
    # is_result_path_b10 = '{}B10%.txt'.format(file_url.format(idd))
    # is_result_path_b20 = '{}B20%.txt'.format(file_url.format(idd))
    # is_result_path_b60 = '{}B60%.txt'.format(file_url.format(idd))
    # is_result_path_b90 = '{}B90%.txt'.format(file_url.format(idd))
    # is_result_path_c10 = '{}C10%.txt'.format(file_url.format(idd))
    # is_result_path_c20 = '{}C20%.txt'.format(file_url.format(idd))
    # is_result_path_c60 = '{}C60%.txt'.format(file_url.format(idd))
    # is_result_path_c90 = '{}C90%.txt'.format(file_url.format(idd))
    # is_result_path_100 = '{}100%.txt'.format(file_url.format(idd))
    # aa1 = '{}aa1.html'.format(file_url.format(idd))
    # aa2 = '{}aa2.html'.format(file_url.format(idd))
    # aa3 = '{}aa3.html'.format(file_url.format(idd))
    # gRNA1 = '{}gRNA1.csv'.format(file_url.format(idd))
    # gRNA2 = '{}gRNA2.csv'.format(file_url.format(idd))
    # gRNA3 = '{}gRNA3.csv'.format(file_url.format(idd))
    # try:
    #     os.remove(is_result_path_a20)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_a60)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_a90)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_b20)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_b60)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_b90)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_c20)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_c60)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_c90)
    # except:
    #     pass
    # try:
    #     os.remove(aa1)
    # except:
    #     pass
    # try:
    #     os.remove(aa2)
    # except:
    #     pass
    # try:
    #     os.remove(aa3)
    # except:
    #     pass
    # try:
    #     os.remove(gRNA1)
    # except:
    #     pass
    # try:
    #     os.remove(gRNA2)
    # except:
    #     pass
    # try:
    #     os.remove(gRNA3)
    # except:
    #     pass
    # try:
    #     os.remove(is_result_path_100)
    # except:
    #     pass





