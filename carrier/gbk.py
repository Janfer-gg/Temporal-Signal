import re
from Bio.Seq import complement


def single_gbk(vector_name,scaffold_name,gRNA):
    with open(r"C:\Users\41518\Desktop\载体\载体图谱\单gRNA\{}.gbk".format(scaffold_name), 'r') as aa:
        gbk_list = []
        for row in aa:
            gbk_list.append(row)
        pos_origin = gbk_list.index('ORIGIN\n')
        origin = gbk_list[pos_origin + 1:-1]
        seq = ''
        for item in origin:
            item = re.sub(r'[^A-Za-z]', "", item)
            seq = seq + item


        name_pos = re.findall("gRNA\d+", vector_name)
        # if not name_pos:
        #     return ''

        gRNA_name = name_pos[0]

        pos_gRNA = gbk_list.index('                     /label=gRNA1\n')
        gbk_list[pos_gRNA]=re.sub("gRNA1",gRNA_name,gbk_list[pos_gRNA])
        pos_gRNA = re.findall(r'\d+', gbk_list[pos_gRNA - 1])

        for item in gbk_list:
            keyword = re.search('KEYWORDS', item)
            if keyword:
                pos_keyword=gbk_list.index(item)
                break

        gbk_list[pos_keyword]='KEYWORDS    '+vector_name +'\n'

        seq = seq[0:int(pos_gRNA[0])] + gRNA + seq[int(pos_gRNA[1]):]

    aa.close()

    gbk_list = gbk_list[0:pos_origin + 1]
    gbk_list.append(seq)
    gbk_list.append("\n//\n")

    f = open("C:/Users/41518/Desktop/{}.dna".format(vector_name), "w", encoding='utf-8')
    f.writelines(gbk_list)
    f.close()


def dual_gbk(vector_name,scaffold_name,gRNA1,gRNA2):
    with open(r"C:\Users\41518\Desktop\载体\载体图谱\双gRNA\{}.gbk".format(scaffold_name), 'r') as aa:
        gbk_list = []
        for row in aa:
            gbk_list.append(row)

        pos_origin = gbk_list.index('ORIGIN\n')
        origin = gbk_list[pos_origin + 1:-1]
        seq = ''
        for item in origin:
            item = re.sub(r'[^A-Za-z]', "", item)
            seq = seq + item

        name_pos = re.findall("gRNA\d+", vector_name)
        gRNA1_name = name_pos[0]
        gRNA2_name = name_pos[1]

        pos_gRNA1 = gbk_list.index('                     /label=gRNA1\n')
        pos_gRNA2=gbk_list.index('                     /label=gRNA2\n')
        gbk_list[pos_gRNA1]=re.sub("gRNA1",gRNA1_name,gbk_list[pos_gRNA1])
        gbk_list[pos_gRNA2] = re.sub("gRNA2", gRNA2_name, gbk_list[pos_gRNA2])

        pos_gRNA1 = re.findall(r'\d+', gbk_list[pos_gRNA1 - 1])
        pos_gRNA2 = re.findall(r'\d+', gbk_list[pos_gRNA2 - 1])

        for item in gbk_list:
            keyword = re.search('KEYWORDS', item)
            if keyword:
                pos_keyword=gbk_list.index(item)
                break

        gbk_list[pos_keyword]='KEYWORDS    '+vector_name +'\n'

        seq = seq[0:int(pos_gRNA1[0])] + gRNA1 + seq[int(pos_gRNA1[1]):int(pos_gRNA2[0])]+ gRNA2 +seq[int(pos_gRNA2[1]):]

    aa.close()

    gbk_list = gbk_list[0:pos_origin + 1]
    gbk_list.append(seq)
    gbk_list.append("\n//\n")

    f = open("C:/Users/41518/Desktop/{}.dna".format(vector_name), "w", encoding='utf-8')
    f.writelines(gbk_list)
    f.close()


def sh_gbk(vector_name,scaffold_name,shRNA):
    with open(r"C:\Users\41518\Desktop\载体\载体图谱\shRNA\{}.gbk".format(scaffold_name), 'r') as aa:
        gbk_list = []
        for row in aa:
            gbk_list.append(row)

        pos_origin = gbk_list.index('ORIGIN\n')
        origin = gbk_list[pos_origin + 1:-1]
        seq = ''
        for item in origin:
            item = re.sub(r'[^A-Za-z]', "", item)
            seq = seq + item

        name_pos = re.findall("shRNA\d+", vector_name)
        gRNA_name = name_pos[0]

        pos_gRNA = gbk_list.index('                     /label=shRNA1\n')
        gbk_list[pos_gRNA]=re.sub("shRNA1",gRNA_name,gbk_list[pos_gRNA])
        pos_gRNA = re.findall(r'\d+', gbk_list[pos_gRNA - 1])

        shRNA_rev = complement(shRNA[::-1])

        for item in gbk_list:
            keyword = re.search('KEYWORDS', item)
            if keyword:
                pos_keyword=gbk_list.index(item)
                break

        gbk_list[pos_keyword]='KEYWORDS    '+vector_name +'\n'

        seq = seq[0:int(pos_gRNA[0])-1] + shRNA+'CTCGAG'+shRNA_rev+'TTTTTG' + seq[int(pos_gRNA[1]):]

    aa.close()

    gbk_list = gbk_list[0:pos_origin + 1]
    gbk_list.append(seq)
    gbk_list.append("\n//\n")

    f = open("C:/Users/41518/Desktop/{}.dna".format(vector_name), "w", encoding='utf-8')
    f.writelines(gbk_list)
    f.close()



# if __name__ == '__main__':
#     gRNA1 = "MMMMMMMMMMMMMMMMMMMM"
#     gRNA2 = 'MMMMMMMMMMMMMMMMMMMM'
#     frame_msg = 'YKO-LV005'
#     plam_name = 'YKO-LV005-NFE2L2[gRNA3-gRNA4]'
#     dual_gbk(plam_name,frame_msg,gRNA1,gRNA2)


if __name__ == '__main__':
    gRNA = "GGGGGGGGGGGGGGGGGGGGG"
    frame_msg = 'YKO-LV011'
    plam_name = 'YKO-LV011-NFE2L2[gRNA1]'
    single_gbk(plam_name,frame_msg,gRNA)