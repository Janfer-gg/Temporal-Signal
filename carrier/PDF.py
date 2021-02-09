import pymysql
import re
from reportlab.lib.enums import TA_RIGHT, TA_LEFT, TA_CENTER
from reportlab.pdfgen.canvas import Canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.units import cm
from reportlab.pdfbase.ttfonts import TTFont

pdfmetrics.registerFont(TTFont('msyh', 'msyh.ttf'))
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle, _baseFontNameB
from reportlab.lib import colors
from reportlab.platypus import Paragraph, Spacer, Image, Table, TableStyle, BaseDocTemplate, \
    PageTemplate, Frame

from svglib.svglib import svg2rlg
from Ubigene.carrier.vectorbuild2 import get_svg
# 页眉
def header(canvas, doc):
    canvas.saveState()

    p = Paragraph("<img src='%s' width='%d' height='%d'/>" % (r'C:\Users\41518\Desktop\work\Ubigene\carrier\head_pdf.png', 560, 28),
                  style=getSampleStyleSheet()['Normal'])  # 使用一个Paragraph Flowable存放图片
    w, h = p.wrap(doc.width, doc.bottomMargin)
    p.drawOn(canvas, doc.leftMargin, 790)
    # canvas.line(doc.leftMargin, doc.bottomMargin + doc.height + 0.5 * cm, doc.leftMargin + doc.width,
    #             doc.bottomMargin + doc.height + 0.5 * cm)  # 画一条横线
    canvas.restoreState()

# 设置页脚
def footer(canvas, doc):

    canvas.saveState()
    p = Paragraph("<img src='%s' width='%d' height='%d'/>" % (r'C:\Users\41518\Desktop\work\Ubigene\carrier\foot_pdf.jpg', 560, 40),
                  style=getSampleStyleSheet()['Normal'])  # 使用一个Paragraph Flowable存放图片
    w, h = p.wrap(doc.width, doc.bottomMargin)
    p.drawOn(canvas, doc.leftMargin, h)
    # 页码
    # pageNumber = canvas.getPageNumber()
    # canvas.drawString(10 * cm, cm, str(pageNumber))
    canvas.restoreState()


def single_rpt(contract,scaffold,name,seq):

    # 基因名
    name_split = re.search(r"{}-(.*)\[(.*)\]".format(scaffold), name, re.I)
    gene_name = name_split.group(1)
    gRNA_label = re.findall('gRNA\d+|shRNA\d+',name)[0]

    # 酶切方案和电泳图谱
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8')
    cursor = db.cursor()
    cursor.execute("SELECT * FROM vector_reporter WHERE type = 'single gRNA'  AND vector = '%s'  " % (scaffold))
    result = cursor.fetchone()
    Digestion = result[4].split(';')
    Electro = result[5]
    picture = Image(r'C:\Users\41518\Desktop\载体\载体报告\{}.png'.format(scaffold))
    picture.drawWidth = 150
    picture.drawHeight = 200

    story = []
    stylesheet = getSampleStyleSheet()
    # 定义各类文本样式
    normalStyle = ParagraphStyle(name='style1',
                                 fontSize=10,
                                 fontName='msyh',
                                 leading=12,
                                 alignment=TA_RIGHT
                                 )

    titleStyle = stylesheet['Title']
    titleStyle.fontSize = 16
    headingStyle = stylesheet['Heading4']
    headingStyle.fontName = _baseFontNameB
    bodyStyle = stylesheet['BodyText']
    bodyStyle.leading=18
    bodyStyle.fontName ='msyh'
    sequence_style = stylesheet['BodyText']
    sequence_style.leading = 16
    #svg图谱
    get_svg('single gRNA', scaffold, name)
    drawing = svg2rlg('C://Users//41518//Desktop//work//Ubigene//carrier//pdf_vector.svg')
    drawing.width, drawing.height = drawing.minWidth() * 0.5, drawing.height * 0.5
    drawing.scale(0.5, 0.5)
    # drawing = Image('C://Users//41518//Desktop//微信图片_20210203105614.png')
    # story.append(drawing)

    # 标题
    content1 = "<para><b><font fontSize=20>gRNA Plasmind Construction Report</font></b> </para>"
    # 协议号
    content2 = "<para><font fontSize=13>Contract ID:</font> <u color='black'><font " \
               "fontSize=13>{}</font></u></para> ".format(contract)
    # 序列
    plasmid_seq = "<para><font backcolor='yellow'>GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC</font>" \
                  "<font color='red'>G{}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC</font> </para>".format(seq)


    # 表格内容
    component_data = [[Paragraph(content1, titleStyle), "", "","", ""],
                      [Paragraph(content2, normalStyle), "","", "", ""],
                      ["Gene", gene_name,"", "gRNA", seq],
                      ["Plasmid", name, "","", ""],
                      [[Paragraph("Plasmid Map", headingStyle), drawing], "","", "", ""],
                      [[Paragraph("Plasmid sequence", headingStyle),
                        Paragraph("Yellow: U6 promoter", bodyStyle),
                        Paragraph("Red: "+gRNA_label, bodyStyle),
                        Paragraph(plasmid_seq,sequence_style)], "", "", "",""],
                      ['Digestion validation','',Digestion[0]+'\n\n'+Digestion[1],'',""],
                      ['Electrophoresis result','',[picture,Paragraph(Electro,bodyStyle)],"",'']
                      ]

    # 创建表格对象，并设定各列宽度
    component_table = Table(component_data, colWidths=[60, 70,60, 60, 230], rowHeights=[40, 30, 30, 30, 300, 180,50,250])
    # 添加表格样式
    component_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (3, -1), 'msyh'),  # 字体
        ('FONTSIZE', (0, 0), (3, -1), 10),  # 字体大小
        ('SPAN', (0, 0), (4, 0)),  # 合并第一行
        ('SPAN', (0, 1), (4, 1)),  # 合并第二行
        # ('BACKGROUND',(0,0),(-1,0), colors.lightskyblue),#设置第一行背景颜色
        ('SPAN', (1, 2), (2, 2)),  # 合并第三行后二三列
        ('SPAN', (1, 3), (4, 3)),  # 合并第四行后四列
        ('SPAN', (0, 4), (4, 4)),  # 合并第五行
        ('SPAN', (0, 5), (4, 5)),  # 合并第六行
        ('SPAN', (0, 6), (1, 6)),  # 合并第七行前两列
        ('SPAN', (2, 6), (4, 6)),  # 合并第七行后三列
        ('SPAN', (0, 7), (1, 7)),  # 合并第八行前两列
        ('SPAN', (2, 7), (4, 7)),  # 合并第八行后三列
        ('ALIGN', (0, 1), (4, 1), 'LEFT'),  # 对齐
        ('ALIGN', (0, 2), (0, 2), 'CENTER'),  # 对齐
        ('ALIGN', (3, 2), (3, 2), 'CENTER'),  # 对齐
        ('ALIGN', (0, 3), (0, 3), 'CENTER'),  # 对齐
        ('ALIGN', (0, 6), (1, 6), 'CENTER'),  # 对齐
        ('ALIGN', (0, 7), (1, 7), 'CENTER'),  # 对齐
        ('VALIGN', (0, 0), (4, 3), 'MIDDLE'),  # 对齐
        ('VALIGN', (0, 4), (4, 5), 'TOP'),  # 对齐
        ('VALIGN', (0, 6), (4, 7), 'MIDDLE'),  # 对齐

        ('LINEBEFORE', (0, 0), (0, 1), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3
        ('LINEAFTER', (0, 0), (-1, 1), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3
        # ('TEXTCOLOR',(0,1),(-2,-1),colors.royalblue),#设置表格内文字颜色
        ('GRID', (0, 2), (-1, -1), 0.3, colors.black),  # 设置表格框线为红色，线宽为0.3
        ('LINEABOVE', (0, 0), (-1, 0), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3

    ]))

    story.append(component_table)

    doc = BaseDocTemplate('单gRNA.pdf', leftMargin=0.6 * cm, rightMargin=0.6 * cm, )

    frame_footer = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
    template = PageTemplate(id='story', frames=frame_footer, onPage=header, onPageEnd=footer,)
    doc.addPageTemplates([template])

    doc.build(story)

def sh_rpt(contract,scaffold,name,seq):

    # 基因名
    name_split = re.search(r"{}-(.*)\[.*\]".format(scaffold), name, re.I)
    gene_name = name_split.group(1)
    gRNA_label = re.findall('gRNA\d+|shRNA\d+',name)[0]

    # 酶切方案和电泳图谱
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8')
    cursor = db.cursor()
    cursor.execute("SELECT * FROM vector_reporter WHERE type = 'shRNA'  AND vector = '%s'  " % (scaffold))
    result = cursor.fetchone()
    Digestion = result[4].split(';')
    Electro = result[5]
    picture = Image(r'C:\Users\41518\Desktop\载体\载体报告\{}.png'.format(scaffold))
    picture.drawWidth = 150
    picture.drawHeight = 200
    story = []
    stylesheet = getSampleStyleSheet()
    # 定义各类文本样式
    normalStyle = ParagraphStyle(name='style1',
                                 fontSize=10,
                                 fontName='msyh',
                                 leading=12,
                                 alignment=TA_RIGHT
                                 )

    titleStyle = stylesheet['Title']
    titleStyle.fontSize = 16
    headingStyle = stylesheet['Heading4']
    headingStyle.fontName = _baseFontNameB
    bodyStyle = stylesheet['BodyText']
    bodyStyle.leading=18
    bodyStyle.fontName ='msyh'
    sequence_style = stylesheet['BodyText']
    sequence_style.leading = 16
    #svg图谱
    get_svg('shRNA', scaffold, name)
    drawing = svg2rlg('C://Users//41518//Desktop//work//Ubigene//carrier//pdf_vector.svg')
    drawing.width, drawing.height = drawing.minWidth() * 0.5, drawing.height * 0.5
    drawing.scale(0.5, 0.5)
    # story.append(drawing)

    # 标题
    content1 = "<para><b><font fontSize=20>gRNA Plasmind Construction Report</font></b> </para>"
    # 协议号
    content2 = "<para><font fontSize=13>Contract ID:</font> <u color='black'><font " \
               "fontSize=13>{}</font></u></para> ".format(contract)

    def DNA_rev_complement(DNA_sequence):
        reverse = DNA_sequence[::-1]
        complement = ''
        for nt in reverse:
            if nt == 'A':
                complement += 'T'
            if nt == 'T':
                complement += 'A'
            if nt == 'G':
                complement += 'C'
            if nt == 'C':
                complement += 'G'
        return complement
    # 序列
    plasmid_seq = "<para><font backcolor='yellow'>GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC</font>" \
                  "<font color='black'>GG</font><font color='red'>{}CTCGAG{}TTTTTG</font> </para>".format(seq,DNA_rev_complement(seq))


    # 表格内容
    component_data = [[Paragraph(content1, titleStyle), "", "","", ""],
                      [Paragraph(content2, normalStyle), "","", "", ""],
                      ["Gene", gene_name,"", "gRNA", seq],
                      ["Plasmid", name, "","", ""],
                      [[Paragraph("Plasmid Map", headingStyle), drawing], "","", "", ""],
                      [[Paragraph("Plasmid sequence", headingStyle),
                        Paragraph("Yellow: U6 promoter", bodyStyle),
                        Paragraph("Red: "+ gRNA_label, bodyStyle),
                        Paragraph(plasmid_seq,sequence_style)], "", "", "",""],
                      ['Digestion validation','',Digestion[0]+'\n\n'+Digestion[1],'',""],
                      ['Electrophoresis result','',[picture,Paragraph(Electro,bodyStyle)],"",'']
                      ]

    # 创建表格对象，并设定各列宽度
    component_table = Table(component_data, colWidths=[60, 70,60, 60, 230], rowHeights=[40, 30, 30, 30, 300, 180,50,250])
    # 添加表格样式
    component_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (3, -1), 'msyh'),  # 字体
        ('FONTSIZE', (0, 0), (3, -1), 10),  # 字体大小
        ('SPAN', (0, 0), (4, 0)),  # 合并第一行
        ('SPAN', (0, 1), (4, 1)),  # 合并第二行
        # ('BACKGROUND',(0,0),(-1,0), colors.lightskyblue),#设置第一行背景颜色
        ('SPAN', (1, 2), (2, 2)),  # 合并第三行后二三列
        ('SPAN', (1, 3), (4, 3)),  # 合并第四行后四列
        ('SPAN', (0, 4), (4, 4)),  # 合并第五行
        ('SPAN', (0, 5), (4, 5)),  # 合并第六行
        ('SPAN', (0, 6), (1, 6)),  # 合并第七行前两列
        ('SPAN', (2, 6), (4, 6)),  # 合并第七行后三列
        ('SPAN', (0, 7), (1, 7)),  # 合并第八行前两列
        ('SPAN', (2, 7), (4, 7)),  # 合并第八行后三列
        ('ALIGN', (0, 1), (4, 1), 'LEFT'),  # 对齐
        ('ALIGN', (0, 2), (0, 2), 'CENTER'),  # 对齐
        ('ALIGN', (3, 2), (3, 2), 'CENTER'),  # 对齐
        ('ALIGN', (0, 3), (0, 3), 'CENTER'),  # 对齐
        ('ALIGN', (0, 6), (1, 6), 'CENTER'),  # 对齐
        ('ALIGN', (0, 7), (1, 7), 'CENTER'),  # 对齐
        ('VALIGN', (0, 0), (4, 3), 'MIDDLE'),  # 对齐
        ('VALIGN', (0, 4), (4, 5), 'TOP'),  # 对齐
        ('VALIGN', (0, 6), (4, 7), 'MIDDLE'),  # 对齐

        ('LINEBEFORE', (0, 0), (0, 1), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3
        ('LINEAFTER', (0, 0), (-1, 1), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3
        # ('TEXTCOLOR',(0,1),(-2,-1),colors.royalblue),#设置表格内文字颜色
        ('GRID', (0, 2), (-1, -1), 0.3, colors.black),  # 设置表格框线为红色，线宽为0.3
        ('LINEABOVE', (0, 0), (-1, 0), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3

    ]))

    story.append(component_table)

    doc = BaseDocTemplate('shRNA.pdf', leftMargin=0.6 * cm, rightMargin=0.6 * cm, )

    frame_footer = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
    template = PageTemplate(id='story', frames=frame_footer, onPage=header, onPageEnd=footer,)
    doc.addPageTemplates([template])

    doc.build(story)


def dual_rpt(contract,scaffold,name,seq1,seq2):

    # 基因名
    name_split = re.search(r"{}-(.*)\[(.*)\]".format(scaffold), name, re.I)
    gene_name = name_split.group(1)
    gRNA_label = re.findall('gRNA\d+|shRNA\d+',name)

    # 酶切方案和电泳图谱
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8')
    cursor = db.cursor()
    cursor.execute("SELECT * FROM vector_reporter WHERE type = 'dual gRNA'  AND vector = '%s'  " % (scaffold))
    result = cursor.fetchone()
    Digestion = result[4].split(';')
    Electro = result[5]
    picture = Image(r'C:\Users\41518\Desktop\载体\载体报告\{}-d.png'.format(scaffold))
    picture.drawWidth = 150
    picture.drawHeight = 200

    story = []
    stylesheet = getSampleStyleSheet()
    # 定义各类文本样式
    normalStyle = ParagraphStyle(name='style1',
                                 fontSize=10,
                                 fontName='msyh',
                                 leading=12,
                                 alignment=TA_RIGHT
                                 )

    titleStyle = stylesheet['Title']
    titleStyle.fontSize = 16
    headingStyle = stylesheet['Heading4']
    headingStyle.fontName = _baseFontNameB
    bodyStyle = stylesheet['BodyText']
    bodyStyle.leading=18
    bodyStyle.fontName ='msyh'
    sequence_style = stylesheet['BodyText']
    sequence_style.leading = 16
    #svg图谱
    get_svg('dual gRNA', scaffold, name)
    drawing = svg2rlg('C://Users//41518//Desktop//work//Ubigene//carrier//pdf_vector.svg')
    drawing.width, drawing.height = drawing.minWidth() * 0.5, drawing.height * 0.5
    drawing.scale(0.5, 0.5)
    # story.append(drawing)

    # 标题
    content1 = "<para><b><font fontSize=20>gRNA Plasmind Construction Report</font></b> </para>"
    # 协议号
    content2 = "<para><font fontSize=13>Contract ID:</font> <u color='black'><font " \
               "fontSize=13>{}</font></u></para> ".format(contract)
    # 序列
    plasmid_seq = "<para><font backcolor='yellow'>GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC</font>" \
                  "<font color='blue'>G{}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC</font><font color='black'>TTTTTTGGCCGGCC</font><font backcolor='yellow'>GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACC</font>" \
                  "<font color='red'>G{}GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC</font> </para>".format(seq1,seq2)


    # 表格内容
    component_data = [[Paragraph(content1, titleStyle), "", "","", ""],
                      [Paragraph(content2, normalStyle), "","", "", ""],
                      ["Gene", gene_name,"", "gRNA", Paragraph('gRNA1：'+seq1+'\n'+'gRNA2：'+seq2,bodyStyle)],
                      ["Plasmid", name, "","", ""],
                      [[Paragraph("Plasmid Map", headingStyle), drawing], "","", "", ""],
                      [[Paragraph("Plasmid sequence", headingStyle),
                        Paragraph("Yellow: U6 promoter", bodyStyle),
                        Paragraph("Blue: " + gRNA_label[0] +"&nbsp&nbsp&nbsp&nbsp Red: " + gRNA_label[1], bodyStyle),
                        Paragraph(plasmid_seq,sequence_style)], "", "", "",""],
                      ['Digestion validation','',Digestion[0]+'\n\n'+Digestion[1],'',""],
                      ['Electrophoresis result','',[picture,Paragraph(Electro,bodyStyle)],"",'']
                      ]

    # 创建表格对象，并设定各列宽度
    component_table = Table(component_data, colWidths=[60, 70, 60, 60, 230], rowHeights=[40, 30, 40, 30, 280, 260, 50, 250])
    # 添加表格样式
    component_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (3, -1), 'msyh'),  # 字体
        ('FONTSIZE', (0, 0), (3, -1), 10),  # 字体大小
        ('SPAN', (0, 0), (4, 0)),  # 合并第一行
        ('SPAN', (0, 1), (4, 1)),  # 合并第二行
        # ('BACKGROUND',(0,0),(-1,0), colors.lightskyblue),#设置第一行背景颜色
        ('SPAN', (1, 2), (2, 2)),  # 合并第三行后二三列
        ('SPAN', (1, 3), (4, 3)),  # 合并第四行后四列
        ('SPAN', (0, 4), (4, 4)),  # 合并第五行
        ('SPAN', (0, 5), (4, 5)),  # 合并第六行
        ('SPAN', (0, 6), (1, 6)),  # 合并第七行前两列
        ('SPAN', (2, 6), (4, 6)),  # 合并第七行后三列
        ('SPAN', (0, 7), (1, 7)),  # 合并第八行前两列
        ('SPAN', (2, 7), (4, 7)),  # 合并第八行后三列
        ('ALIGN', (0, 1), (4, 1), 'LEFT'),  # 对齐
        ('ALIGN', (0, 2), (0, 2), 'CENTER'),  # 对齐
        ('ALIGN', (3, 2), (3, 2), 'CENTER'),  # 对齐
        ('ALIGN', (0, 3), (0, 3), 'CENTER'),  # 对齐
        ('ALIGN', (0, 6), (1, 6), 'CENTER'),  # 对齐
        ('ALIGN', (0, 7), (1, 7), 'CENTER'),  # 对齐
        ('VALIGN', (0, 0), (4, 3), 'MIDDLE'),  # 对齐
        ('VALIGN', (0, 4), (4, 5), 'TOP'),  # 对齐
        ('VALIGN', (0, 6), (4, 7), 'MIDDLE'),  # 对齐

        ('LINEBEFORE', (0, 0), (0, 1), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3
        ('LINEAFTER', (0, 0), (-1, 1), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3
        # ('TEXTCOLOR',(0,1),(-2,-1),colors.royalblue),#设置表格内文字颜色
        ('GRID', (0, 2), (-1, -1), 0.3, colors.black),  # 设置表格框线为红色，线宽为0.3
        ('LINEABOVE', (0, 0), (-1, 0), 0.3, colors.black),  # 设置表格左边线颜色为灰色，线宽为0.3

    ]))

    story.append(component_table)

    doc = BaseDocTemplate('双gRNA.pdf', leftMargin=0.6 * cm, rightMargin=0.6 * cm, )

    frame_footer = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
    template = PageTemplate(id='story', frames=frame_footer, onPage=header, onPageEnd=footer,)
    doc.addPageTemplates([template])

    doc.build(story)



if __name__ == '__main__':

    # contract = 'UBISVK201028HJ1'
    # scaffold = 'YKO-LV001'
    # name = 'YKO-LV001-BRAF(c.T1799A)[gRNA1-gRNA2]'
    # seq = 'GCATCCTACTTGTCCTAATAA'
    # seq2 = 'TAATGCAGAGTGTGCTCGTC'
    # dual_rpt(contract, scaffold, name, seq, seq2)

    contract = 'UBISVK201028HJ1'
    scaffold = 'YKO-LV001'
    name = 'YKO-LV001-BRAF(c.T1799A)[gRNA3]'
    seq = 'GCATCCTACTTGTCCTAATAA'
    single_rpt(contract, scaffold, name, seq)

    # contract = 'UBISVK201028HJ1'
    # scaffold = 'YSH-LV001'
    # name = 'YSH-LV001-BRAF(c.T1799A)[shRNA5]'
    # seq = 'GCATCCTACTTGTCCTAATAA'
    # sh_rpt(contract, scaffold, name, seq)


