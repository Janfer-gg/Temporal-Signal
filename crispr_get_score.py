import requests
import json
from lxml import etree
import time
import csv
from requests.adapters import HTTPAdapter

requests.packages.urllib3.disable_warnings()

s = requests.Session()
s.mount('http://', HTTPAdapter(max_retries=100))
s.mount('https://', HTTPAdapter(max_retries=100))
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.97 Safari/537.36",
    "Host":"cctop.cos.uni-heidelberg.de",
    "Origin":"https://cctop.cos.uni-heidelberg.de",
    "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3",
    "Referer":"https://cctop.cos.uni-heidelberg.de/",
}


def get_redirect_url(gene_list,species):

    species_dict = {
        'Human':"hg38",
        'Mouse':"mm10",
        'Rat':"rattus_norvegicus_V99",
        'dog':"canis_familiaris_V99",
        'chinese Hamster':"cricetulus_griseus_crigri_V99",
        'green monkey':'green monkey'


    }
    FORM_DATA = {
        "name": 'unnamed',
        "radQ": 'single',
        "sequence": gene_list,
        "pamType": 'NGG',
        "targetLength": '20',
        "sgRNA5": 'NN',
        "sgRNA3": 'NN',
        "inVitroTx": 'T7',
        "totalMismatches": '4',
        "useCore": 'on',
        "coreLength": '12',
        "coreMismatches": '2',
        "species": species_dict[species],
    }
    URL = 'https://cctop.cos.uni-heidelberg.de/cgi-bin/search.py'
    response = s.post(URL, headers=headers, data=FORM_DATA, allow_redirects=False, verify=False)
    new_url = response.headers["Location"]
    return new_url


def get_content(SID):
    THE_RESULT_URL = f'https://cctop.cos.uni-heidelberg.de/result/{SID}/unnamed_frame.html'
    # 执行10次每次停留3s
    for i in range(30):
        time.sleep(3)
        RESULE_RESP = s.get(url=THE_RESULT_URL, verify=False, headers=headers)
        content = RESULE_RESP.text
        score = etree.HTML(content).xpath("//table/tr[2]/td[2]/text()")

        if score:
            score = score[0].split(' ')[0]
            return content


def reader_writer(file_path,species):
    f_r = open(file_path, 'r')
    reader = csv.reader(f_r)
    first_line = next(reader)
    first_line.append('crispr_score')
    f_w = open(file_path, 'w', newline='')
    writer = csv.writer(f_w)
    writer.writerow(first_line)
    read_lines = list(reader)

    gene_list_list = []
    for row in read_lines:
        gene_list = row[2]

        gene_list_list.append(gene_list)

    print(gene_list_list)

    chunk_gene_list = chunkIt(gene_list_list,20)

    Sequence_socre_dict = {}

    for chunk_gene_l in chunk_gene_list:

        chunk_gene_str = ''.join(chunk_gene_l)

        REDIRECT_URL = get_redirect_url(chunk_gene_str,species)
        SID = REDIRECT_URL.replace('search.py?sid=', '')
        content = get_content(SID)

        import re
        Sequence = re.findall(r'<td>Sequence:</td><td class="mono">(.*?)</td>', content)
        socre = re.findall(r'<td bgcolor="(.*?)">(.*?)</td>', content)

        for seq,sco in zip(Sequence,socre):

            Sequence_socre_dict[seq] = sco

    for row in read_lines:
        gg = row[2]

        score = Sequence_socre_dict[gg]

        score = score[1].split('  ')[0].strip()

        print(f"[+]序列{gene_list}的sid={SID}它的 score是{score}")
        row.append(score)
        writer.writerow(row)

    f_r.close()
    f_w.close()


def chunkIt(seq, avg):
    avg = avg
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out
# if __name__ == '__main__':
# #     #TODO 改为自己文件路径
#     file_path = r'C:\Users\41518\Desktop\work\Ubigene\CCTOP-predictor.csv'
#     reader_writer(file_path)
