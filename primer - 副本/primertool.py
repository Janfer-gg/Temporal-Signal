import requests
import time
import re
from lxml import etree
import logging
import csv
import os
import pathlib
from requests.adapters import HTTPAdapter

s = requests.Session()
s.mount('http://', HTTPAdapter(max_retries=100))
s.mount('https://', HTTPAdapter(max_retries=100))

requests.packages.urllib3.disable_warnings()

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)s %(levelname)s %(message)s",
                    datefmt='%Y-%m-%d  %H:%M:%S %a')

result = dict()
gt3_break = list()


class PrimerTool(object):
    def __init__(self, PCR1, PCR3, ORGANISM, MAX_TARGET_SIZE):
        self.headers = {
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3",
            "Host": "www.ncbi.nlm.nih.gov",
            "Origin": "https://www.ncbi.nlm.nih.gov",
            "Referer": "https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?LINK_LOC=BlastHome",
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/78.0.3904.97 Safari/537.36",
        }

        ORGANISM_dict = {
            'Human': 'Homo sapiens',
            'Mouse': 'mouse (taxid:10088)',
            'Rat': 'Norway rat (taxid:10116)',
            'dog': 'dog (taxid:9615)',
            'chinese Hamster': '',
            'Zebrafish': 'zebrafish (taxid:7955)',
            'Chicken': '',
            'Pig': '',
            'Goat': '',
        }
        self.data = {
            "INPUT_SEQUENCE": "",
            "SEQFILE": "",
            "PRIMER5_START": "",
            "PRIMER5_END": "",
            "PRIMER3_START": "",
            "PRIMER3_END": "",
            "PRIMER_LEFT_INPUT": PCR1,
            "PRIMER_RIGHT_INPUT": PCR3,
            "PRIMER_PRODUCT_MIN": "70",
            "PRIMER_PRODUCT_MAX": "1000",
            "PRIMER_NUM_RETURN": "10",
            "PRIMER_MIN_TM": "57.0",
            "PRIMER_OPT_TM": "60.0",
            "PRIMER_MAX_TM": "63.0",
            "PRIMER_MAX_DIFF_TM": "3",
            "PRIMER_ON_SPLICE_SITE": "0",
            "SPLICE_SITE_OVERLAP_5END": "7",
            "SPLICE_SITE_OVERLAP_3END": "4",
            "SPLICE_SITE_OVERLAP_3END_MAX": "8",
            "MIN_INTRON_SIZE": "100",
            "MAX_INTRON_SIZE": "1000000",
            "SEARCH_SPECIFIC_PRIMER": "on",
            "SEARCHMODE": "0",
            "PRIMER_SPECIFICITY_DATABASE": "PRIMERDB/genome_selected_species",
            "CUSTOM_DB": "",
            "CUSTOMSEQFILE": "",
            "ORGANISM": ORGANISM_dict[ORGANISM],
            "ORG_DBS": "on",
            "slctOrg": "",
            "ENTREZ_QUERY": "",
            "TOTAL_PRIMER_SPECIFICITY_MISMATCH": "1",
            "PRIMER_3END_SPECIFICITY_MISMATCH": "1",
            "MISMATCH_REGION_LENGTH": "5",
            "TOTAL_MISMATCH_IGNORE": "6",
            "MAX_TARGET_SIZE": str(MAX_TARGET_SIZE),
            "HITSIZE": "50000",
            "UNGAPPED_BLAST": "on",
            "EVALUE": "30000",
            "WORD_SIZE": "7",
            "MAX_CANDIDATE_PRIMER": "500",
            "NUM_TARGETS": "20",
            "NUM_TARGETS_WITH_PRIMERS": "1000",
            "MAX_TARGET_PER_TEMPLATE": "100",
            "PRODUCT_MIN_TM": "",
            "PRODUCT_OPT_TM": "",
            "PRODUCT_MAX_TM": "",
            "PRIMER_MIN_SIZE": "15",
            "PRIMER_OPT_SIZE": "20",
            "PRIMER_MAX_SIZE": "25",
            "PRIMER_MIN_GC": "20.0",
            "PRIMER_MAX_GC": "80.0",
            "GC_CLAMP": "0",
            "POLYX": "5",
            "PRIMER_MAX_END_STABILITY": "9",
            "PRIMER_MAX_END_GC": "5",
            "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": "40.00",
            "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": "70.00",
            "PRIMER_MAX_SELF_ANY_TH": "45.0",
            "PRIMER_MAX_SELF_END_TH": "35.0",
            "PRIMER_PAIR_MAX_COMPL_ANY_TH": "45.0",
            "PRIMER_PAIR_MAX_COMPL_END_TH": "35.0",
            "PRIMER_MAX_HAIRPIN_TH": "24.0",
            "PRIMER_MAX_TEMPLATE_MISPRIMING": "12.00",
            "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": "24.00",
            "SELF_ANY": "8.00",
            "SELF_END": "3.00",
            "PRIMER_PAIR_MAX_COMPL_ANY": "8.00",
            "PRIMER_PAIR_MAX_COMPL_END": "3.00",
            "EXCLUDED_REGIONS": "",
            "OVERLAP": "",
            "OVERLAP_5END": "7",
            "OVERLAP_3END": "4",
            "MONO_CATIONS": "50.0",
            "DIVA_CATIONS": "1.5",
            "CON_DNTPS": "0.6",
            "SALT_FORMULAR": "1",
            "TM_METHOD": "1",
            "CON_ANEAL_OLIGO": "50.0",
            "PRIMER_MISPRIMING_LIBRARY": "AUTO",
            "LOW_COMPLEXITY_FILTER": "on",
            "PRIMER_INTERNAL_OLIGO_MIN_SIZE": "18",
            "PRIMER_INTERNAL_OLIGO_OPT_SIZE": "20",
            "PRIMER_INTERNAL_OLIGO_MAX_SIZE": "27",
            "PRIMER_INTERNAL_OLIGO_MIN_TM": "57.0",
            "PRIMER_INTERNAL_OLIGO_OPT_TM": "60.0",
            "PRIMER_INTERNAL_OLIGO_MAX_TM": "63.0",
            "PRIMER_INTERNAL_OLIGO_MIN_GC": "20.0",
            "PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT": "50",
            "PRIMER_INTERNAL_OLIGO_MAX_GC": "80.0",
            "LINK_LOC": "BlastHome",
            "SVIEWER_DATA_KEY": '',
            "CMD": "request",
            "NUM_DIFFS": "3",
            "NUM_OPTS_DIFFS": "0",

        }

    def get_job_url(self, **kwargs):
        resp = requests.post(**kwargs)
        headers = resp.headers
        NCBI_RCGI_RetryURL = headers['NCBI-RCGI-RetryURL']
        return NCBI_RCGI_RetryURL

    def get_result(self, NCBI_RCGI_RetryURL):
        # 执行10次每次停留3s
        for i in range(30):
            RESULE_RESP = s.get(url=NCBI_RCGI_RetryURL, verify=False)
            time.sleep(3)
            content = RESULE_RESP.text
            if 'Primer-BLAST Results' in content:
                product_sorce = re.findall(r'product length = (\d+)', content, re.S)

                return product_sorce
                # break


def run_primertool(PCR1_list, PCR3_list, ORGANISM, MAX_TARGET_SIZE):
    start_url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'

    result_list = list()

    for PCR1_row in PCR1_list:
        # print(PCR1)
        for PCR3_row in PCR3_list:
            PCR1, PCR3 = PCR1_row, PCR3_row
            print(PCR1, PCR3)
            primertool = PrimerTool(PCR1=PCR1, PCR3=PCR3, ORGANISM=ORGANISM, MAX_TARGET_SIZE=MAX_TARGET_SIZE)
            NCBI_RCGI_RetryURL = primertool.get_job_url(url=start_url, headers=primertool.headers, data=primertool.data,
                                                        verify=False)
            result_content_list = primertool.get_result(NCBI_RCGI_RetryURL)
            result_content_num = len(result_content_list)
            print(result_content_num)
            result['{},{}'.format(PCR1, PCR3)] = result_content_list
            if int(result_content_num) == 1:
                # result_list.append([PCR1,PCR3])
                return [PCR1, PCR3]

    for reskey, resval in result.items():

        result_content_n = len(resval)

        if int(result_content_n) == 2 and (resval[1] == resval[0]):
            PCR1, PCR3 = reskey.split(',')

            result_list.append([PCR1, PCR3])

            if len(gt3_break) >= 3:
                return result_list[-1]
        else:
            pass

#
# if __name__ == '__main__':
#     # TODO: 更改文件路径
#     PCR1_list = ['GGTGATGCTTTGTTGTTGGTAGC']
#     PCR3_list = ['GGTTGGGACAGCCCAAGTTG', 'CACAGAGGAATGCAGCCCAC', 'GACACAGAGGAATGCAGCCC', 'GACACAGAGGAATGCAGCCCA',
#                  'GGACACAGAGGAATGCAGCC']
#     ORGANISM = 'Human'
#     result = run_primertool(PCR1_list, PCR3_list, ORGANISM, 4000)
#
#     print(result)
