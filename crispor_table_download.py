import requests
from lxml import etree
import re
import logging
import csv
import time

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)s %(levelname)s %(message)s",
                    datefmt='%Y-%m-%d  %H:%M:%S %a')
url = 'http://crispor.tefor.net/crispor.py'
headers = {
    "Host": "crispor.tefor.net",
    "Connection": "keep-alive",
    "Upgrade-Insecure-Requests": "1",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.97 Safari/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3",
    "Referer": "http://crispor.tefor.net/crispor.py",
    "Accept-Language": "zh-CN,zh;q=0.9,en;q=0.8",
}


def post_msg(**kwargs):

    resp = requests.post(**kwargs)
    return resp.text


def get_msg(**kwargs):
    resp = requests.post(**kwargs)
    return resp


def xpath_html(content):
    return etree.HTML(content)


def run(array,org):

    species_dict = {
        'Human':"hg38",
        'Mouse':"mm10",
        'Rat':"ensRatNor",
        'dog':"ens79CanFam",
        'chinese Hamster':"GCA_900186095.1",
        'green monkey':'green monkey'


    }
    data = {
        "name": "",
        "seq": array,
        "org": species_dict[org],
        "pam": "NGG",
        "submit": "SUBMIT",
    }
    content = post_msg(url=url, data=data, headers=headers)
    time.sleep(1)
    #logging.info(content)
    batch_id = re.findall(r"history.replaceState\('crispor.py', document.title, '\?batchId=(.*?)'\);", content)[0]
    logging.info("*" * 100 + batch_id)
    batch_302_url = 'http://crispor.tefor.net/crispor.py?batchId={}'.format(batch_id)
    html_xpath = ''
    while True:
        logging.info(batch_302_url)
        resp = get_msg(url=batch_302_url, headers=headers, stream=True)
        logging.info('This page will refresh every 10 seconds')
        time.sleep(10)
        logging.info('等待结束')
        logging.info('数据正在加载')
        with open('aa.html', 'wb')as f:
            for chunk in resp.iter_content(chunk_size=512):
                if chunk:
                    f.write(chunk)
        logging.info('数据加载完毕')
        with open('aa.html', 'rb')as r_f:
            resp_content = r_f.read()
        html_xpath = xpath_html(resp_content)
        thead_tr_list = html_xpath.xpath("//*[@id='otTable']/thead/tr")
        # trs = html_xpath.xpath("//*//tr[contains(@id,'s')]")
        logging.info(thead_tr_list)
        if thead_tr_list:
            break
    with open('aa.html', 'rb')as r_f:
        resp_content = r_f.read()
    html_xpath = xpath_html(resp_content)

    thead_tr_list = html_xpath.xpath("//*[@id='otTable']/thead/tr")
    trs = html_xpath.xpath("//tr")

    with open('gRNA.csv'.format(batch_id), 'w', encoding='utf8', newline='')as aa:
        csv_writer = csv.writer(aa)
        for thead_tr in thead_tr_list:
            thead_ths = thead_tr.xpath('./th')
            thead_ths_list = list()
            for thead_td in thead_ths:
                th_str = thead_td.xpath('string(.)').strip()
                thead_ths_list.append(th_str)
            print(thead_ths_list)
        for tr in trs:
            tds = tr.xpath('./td')
            td_list = list()
            for td in tds:
                td_str = td.xpath('string(.)').strip()
                td_list.append(td_str)
            # logging.info('[+]正在写入行{}'.format(td_list))
            csv_writer.writerow(td_list)
        logging.debug('[+]写入gRNA.csv完毕'.format(batch_id))

#
# if __name__ == '__main__':
#     array = 'ATTTTCTCAAGTTGACATTAGAGAACTACGTGAAAACTTTAGAAGTGCCAGTGAAGATTATCAGGATGGAAGAGCGCTCTGGGTTAATACGCGCCCGTCTCCGAGGAGCAGCAGCTTCAAAAGGACAGGTCATAACTTTTCTGGATGCACACTGCGAGTGCACCTTAGGATGGCTGGAGCCCTTGCTGGCGAGAATAAAGGAAGACAG'
#     run(array,'Mouse')
