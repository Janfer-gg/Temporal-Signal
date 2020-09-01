import requests
import json
from lxml import etree
import csv
import logging
import re
logging.basicConfig(level=logging.DEBUG,format="%(asctime)s %(name)s %(levelname)s %(message)s",datefmt = '%Y-%m-%d  %H:%M:%S %a' )
headers = {
"Accept":"application/json, text/javascript, */*; q=0.01",
"Accept-Encoding":"gzip, deflate",
"Accept-Language":"zh-CN,zh;q=0.9,en;q=0.8",
"Connection":"keep-alive",
"Cookie":"redirect_mirror=no; ENSEMBL_WWW_SESSION=056b9651e8fb71e20f0121c65f49f0fe0b1c499d3349697c7d2083e4; DYNAMIC_WIDTH=1; _ga=GA1.2.493265046.1593497440; _gid=GA1.2.316210055.1593497440; ENSEMBL_WIDTH=1600; _gat=1; _gali=zmenu_150_5_1596_457",
"Host":"asia.ensembl.org",
"Referer":"http://asia.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000124523;r=6:13574529-13615158",
"User-Agent":"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.97 Safari/537.36",
"X-Requested-With":"XMLHttpRequest",
}
def get_msg(**kwargs):
    resp = requests.get(**kwargs)
    return resp.text
def xpath_html(content):
    return etree.HTML(content)

def download_csv(gene_name):
    start_url = 'http://asia.ensembl.org/Gene/Summary?db=core;g={gene_name}'.format(gene_name=gene_name)
    resp_content = get_msg(url=start_url,headers=headers)
    resp_content = re.sub(r'<span class="_ht_tip hidden">(.*?)</span>','',resp_content)
    html_xpath = xpath_html(resp_content)
    th_list = html_xpath.xpath('//*[@id="transcripts_table"]//th//text()')
    trs = html_xpath.xpath('//*[@id="transcripts_table"]//tbody/tr')
    with open('transcript.csv', 'w', encoding='utf8',newline='')as aa:
        csv_writer = csv.writer(aa)
        csv_writer.writerow(th_list)
        for tr in trs:
            tds = tr.xpath('./td')
            td_list = list()
            for td in tds:
                td_str = td.xpath('string(.)').strip()
                td_list.append(td_str)
            td_list[-1]=td_list[-1].split('\n        \n        ')[-1]
            logging.info('[+]正在写入行{}'.format(td_list))
            csv_writer.writerow(td_list)
        logging.debug('[+]写入{}.csv完毕'.format(gene_name))
# if __name__ == '__main__':
#     #TODO 更改gene_name
#     gene_name = 'ENSG00000124523'
#     download_csv(gene_name)
