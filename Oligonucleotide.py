from selenium import webdriver
from fake_useragent import UserAgent
from lxml import etree
import csv
from selenium.webdriver.support.wait import WebDriverWait
from  selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import time
import re
from selenium.webdriver.common.keys import Keys


class AddClickNum(object):


    def __init__(self,lujing,google_zhushou):

        ua = UserAgent()
        self.time = 2
        self.searh_page_time = 12
        chromeOptions = webdriver.ChromeOptions()
        # chromeOptions.add_argument('lang=zh_CN.UTF-8')
        chromeOptions.add_argument('--headless')
        # chromeOptions.add_argument('user-agent="{}"'.format(ua))
        # 设置好应用扩展
        extension_path = google_zhushou
        # chromeOptions.add_extension(extension_path)
        self.lujing = lujing
        self.chromeOptions = chromeOptions
    def init_driver(self):
        self.driver = webdriver.Chrome(executable_path=self.lujing,chrome_options=self.chromeOptions)
        self.weizhi = 0
    def search_key(self,keyword):
        self.init_driver()

        import random
        self.driver.get('http://biotools.nubic.northwestern.edu/OligoCalc.html')

        self.driver.find_element_by_name('oligoBox').send_keys(keyword)
        time.sleep(1)
        self.driver.find_element_by_name('Calbutton').click()

        Complement = self.driver.find_element_by_xpath('//*[@onclick="return calcPrimer(this.form)"]')
        Complement.click()


        # WebDriverWait(self.driver, 10).until(EC.alert_is_present())
        # self.driver.switch_to.alert.send_keys('genzhengmiaobuhong@gmail.com')
        # self.driver.switch_to.alert.accept()
        first_handle = self.driver.current_window_handle

        handles = self.driver.window_handles
        for i in handles:
            if i != first_handle:
                self.driver.close()  # 关闭当前窗口
                self.driver.switch_to.window(i)
                page_source = self.driver.page_source
                None_count = re.findall(r'None',page_source)
                len_count = len(None_count)
                print(len_count)



    def quit(self):
        # time.sleep(3600)
        self.driver.quit()
if __name__ == '__main__':

    # todo 改成自己的chromdriver路径
    chrome_driver = r'C:\Program Files (x86)\Google\Chrome\Application\chrome'
    # todo 这里是谷歌访问助手的路径
    google_zhushou = r'E:\谷歌访问助手_v2.3.0.crx'


    add_click_num = AddClickNum(chrome_driver,google_zhushou)
    # file_path = 'E:\workspace\workspace\project1\PCR1.csv'
    # f_r = open(file_path, 'r')
    # reader = csv.reader(f_r)
    # first_line = next(reader)
    # first_line.append('crispr_score')
    # f_w = open(file_path, 'w', newline='')
    # writer = csv.writer(f_w)
    # writer.writerow(first_line)
    # for row in reader:
    #     print(row)
        # keyword = row[0]
        # writer.writerow(row)

    s = time.time()

    for i in ['ATG CAT GCG CTT AGC GTC TA','ATG CAT GCG CTT AGC GTC T','ATG CAT GCG CTT AGC GTC']:
        add_click_num.search_key(keyword=i)
        add_click_num.quit()

    e = time.time()

    t = e - s
    print(t)








