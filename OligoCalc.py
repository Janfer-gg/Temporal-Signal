from selenium import webdriver
from fake_useragent import UserAgent
import csv
import time
import re


class AddClickNum(object):
    def __init__(self, lujing):

        ua = UserAgent()
        self.time = 2
        self.searh_page_time = 12
        chromeOptions = webdriver.ChromeOptions()
        chromeOptions.add_argument('--headless')
        # chromeOptions.add_argument('user-agent="{}"'.format(ua))
        # 设置应用扩展
        # chromeOptions.add_extension(extension_path)
        chrome_prefs = {}

        chromeOptions.experimental_options["prefs"] = chrome_prefs
        chrome_prefs["profile.default_content_settings"] = {"images": 2}
        chrome_prefs["profile.managed_default_content_settings"] = {"images": 2}
        self.lujing = lujing
        self.chromeOptions = chromeOptions
        # def init_driver(self):
        self.driver = webdriver.Chrome(executable_path=self.lujing, chrome_options=self.chromeOptions)
        self.driver.get('http://biotools.nubic.northwestern.edu/OligoCalc.html')
        self.weizhi = 0

    def search_key(self, keyword):
        # self.init_driver()
        import random
        self.driver.get('http://biotools.nubic.northwestern.edu/OligoCalc.html')

        self.driver.find_element_by_name('oligoBox').send_keys(keyword)
        time.sleep(1)
        self.driver.find_element_by_name('Calbutton').click()

        Complement = self.driver.find_element_by_xpath('//*[@onclick="return calcPrimer(this.form)"]')
        Complement.click()

        first_handle = self.driver.current_window_handle

        handles = self.driver.window_handles
        print(handles)
        for i in handles:
            if i != first_handle:
                self.driver.switch_to.window(i)
                page_source = self.driver.page_source
                None_count = re.findall(r'None', page_source)
                len_count = len(None_count)
                self.driver.close()  # 关闭当前窗口
                self.driver.switch_to.window(first_handle)

                return len_count

    def quit(self):
        # time.sleep(3600)
        self.driver.quit()


def RunOligoCalc(file_path, chrome_driver):
    f_r = open(file_path, 'r')
    reader = csv.reader(f_r)
    first_line = next(reader)
    print(first_line)
    print('*' * 100)

    reader_list = list(reader)
    print(reader_list)

    f_w = open(file_path, 'w', newline='')
    writer = csv.writer(f_w)
    writer.writerow(first_line)
    f_w.close()


    add_click_num = AddClickNum(chrome_driver)
    for row in reader_list:

        print(row)
        keyword = row[0]
        none_count = add_click_num.search_key(keyword=keyword)
        print(none_count)
        if int(none_count) == 3:
            f_ww = open(file_path, 'a', newline='')
            wwriter = csv.writer(f_ww)
            wwriter.writerow(row)
            f_ww.close()

    f_r.close()
    add_click_num.quit()


# if __name__ == '__main__':
#     # todo 改成自己的chromdriver路径
#     chrome_driver = r'E:\chromedriver.exe'
#     # todo 改成自己的文件路径
#
#     file_path = 'E:\workspace\workspace\project1\PCR1.csv'
#
#     s = time.time()
#     RunOligoCalc(file_path, chrome_driver)
#     e = time.time()
#     t = e - s
#     print(t)
