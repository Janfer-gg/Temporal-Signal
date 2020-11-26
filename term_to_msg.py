#!/usr/bin/env python
#-*- coding:utf-8 -*-
# author:NAIVE
# datetime:2020/8/21 11:44
import requests
import json
import re
from scrapy.selector import Selector

import logging
from requests.adapters import HTTPAdapter

requests.packages.urllib3.disable_warnings()
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)s %(levelname)s %(message)s",
                    datefmt='%Y-%m-%d  %H:%M:%S %a')
url = 'http://crispor.tefor.net/crispor.py'


s = requests.Session()

s.mount('http://', HTTPAdapter(max_retries=100))
s.mount('https://', HTTPAdapter(max_retries=100))
term_get_id_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}&sort=Relevance&retmode=xml'
id_get_summary_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={}&retmode=xml'
id_get_esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&retmode=xml'


def get_resp_text(url, **kwargs):
    resp = s.get(url=url)
    return resp.text


def xml_selector(api_url):
    xml = get_resp_text(api_url)
    sel = Selector(type='xml', text=xml)

    return sel


def get_id(term):
    sel = xml_selector(api_url=term_get_id_url.format(term))
    id_list = sel.xpath('//Id//text()').extract()

    return id_list


def get_summary(term, scitific):
    id_list = get_id(term)

    Scientific_list = list()
    result_list = list()
    if not scitific:
        for idd in id_list:
            print(idd)

            item, ScientificName = sel_get_summary(idd)

            print(ScientificName)

            Organism_list = ['Cricetulus griseus', 'Homo sapiens', 'Mus musculus', 'Rattus rattus', 'Rattus norvegicus',
                             'Canis lupus familiaris', 'Chlorocebus sabaeus']
            if (ScientificName not in Scientific_list) and (ScientificName in Organism_list):
                Scientific_list.append(ScientificName)
                result_list.append(item)
            else:
                pass
            if len(Scientific_list) > 3:
                break


    else:
        for idd in id_list:
            if len(Scientific_list) > 3:
                break
            item, ScientificName = sel_get_summary(idd)

            Scientific_list.append(idd)

            # item, ScientificName = sel_get_summary(idd)
            result_list.append(item)

    return result_list


def sel_get_summary(idd):
    item = dict()
    sel_summary = xml_selector(api_url=id_get_summary_url.format(idd))
    item['OfficalSymbol'] = sel_summary.xpath('//NomenclatureSymbol//text()').extract_first()

    if not item['OfficalSymbol']:
        item['OfficalSymbol'] = sel_summary.xpath('//Name//text()').extract_first()
    item['OfficalFullName'] = sel_summary.xpath('//NomenclatureName//text()').extract_first()
    item['Organism'] = sel_summary.xpath('//ScientificName//text()').extract_first()
    item['OtherAliases'] = sel_summary.xpath('//OtherAliases//text()').extract_first()
    item['Chromosome'] = sel_summary.xpath('//Chromosome//text()').extract_first()
    try:
        item['Summary'] = re.sub('\[.*?\]', '', sel_summary.xpath('//Summary//text()').extract_first())
    except:
        item['Summary'] = sel_summary.xpath('//Summary//text()').extract_first()

    item['NcbiGeneId'] = idd

    sel_esearch = xml_selector(api_url=id_get_esearch_url.format(idd))

    # esearch_text = get_resp_text(url=id_get_esearch_url.format(idd))
    # esearch_res = re.search(r'type (.*?),',esearch_text).group()
    # item['GeneType'] = esearch_res.replace('type ','').replace(',','')
    EnsemblIdList = sel_esearch.xpath('//Entrezgene_gene//Object-id_str//text()').extract()

    print(EnsemblIdList)
    for i in EnsemblIdList:

        if i.startswith('EN') and len(i) > 9:
            item['EnsemblId'] = i
            break
        else:
            item['EnsemblId'] = ''
    item['Entrezgene_type'] = sel_esearch.xpath('//Entrezgene_type /@value').extract_first()

    print(item)
    return item, item['Organism']


if __name__ == '__main__':
    # import time
    #
    # s = time.time()
    # item_list = get_summary('(PPIA) AND "Homo sapiens"[porgn:__txid9606]',"Homo sapiens")
    # e = time.time()
    # print(e - s)
    # print(len(item_list))
    # print(item_list)
    id_list = get_id('b0001')
    print(id_list)


    # k,v = sel_get_summary(54880)

    # print(k)

    # sel_get_summary(5478)
