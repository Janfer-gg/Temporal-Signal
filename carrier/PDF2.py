from reportlab.pdfgen.canvas import Canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.enums import TA_RIGHT, TA_LEFT, TA_CENTER
from reportlab.lib.units import cm
# from reportlab.pdfbase.cidfonts import UnicodeCIDFont
# pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
from reportlab.pdfbase.ttfonts import TTFont
pdfmetrics.registerFont(TTFont('msyh', 'msyh.ttf'))
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, BaseDocTemplate,PageTemplate, Frame
from reportlab.lib.units import inch
import os
from svglib.svglib import svg2rlg
result_dict=[{'name': 'RSV promoter', 'type': 'Promoter', 'description': 'Rous sarcoma virus enhancer/promoter', 'notes': 'Strong promoter; drives transcription of viral RNA in packaging cells.', 'start': 2, 'end': 229, 'color': '#8F4987'}, {'name': "5' LTR (truncated)", 'type': 'LTR', 'description': "Truncated HIV-1 5' long terminal repeat", 'notes': 'Allows transcription of viral RNA and its packaging into virus.', 'start': 229, 'end': 410, 'color': '#C64394'}, {'name': 'HIV-1 Ψ', 'type': 'Miscellaneous', 'description': 'HIV-1 packaging signal', 'notes': 'Allows packaging of viral RNA into virus.', 'start': 456, 'end': 582, 'color': '#441132'}, {'name': 'RRE', 'type': 'Miscellaneous', 'description': 'HIV-1 Rev response element', 'notes': 'Rev protein binding site that allows Rev-dependent nuclear export of viral RNA during viral packaging.', 'start': 1074, 'end': 1308, 'color': '#237FBB'}, {'name': 'cPPT/CTS', 'type': 'Miscellaneous', 'description': 'Central polypurine tract', 'notes': 'Facilitates the nuclear import of HIV-1 cDNA through a central DNA flap.', 'start': 1802, 'end': 1920, 'color': '#6C9E4F'}, {'name': 'U6 promoter', 'type': 'Promoter', 'description': 'Human U6 small nuclear 1 promoter', 'notes': 'Pol III promoter; drives expression of small RNAs.', 'start': 1926, 'end': 2175, 'color': '#31493A'}, {'name': 'sample', 'type': 'gRNA', 'description': '', 'notes': 'Target: None in human and mouse.', 'start': 2176, 'end': 2270, 'color': '#4A9CEA'}, {'name': 'hPGK', 'type': 'Promoter', 'description': '', 'notes': '', 'start': 2282, 'end': 2787, 'color': '#8F56E7'}, {'name': 'EGFP:T2A:Hygro', 'type': 'Marker', 'description': '', 'notes': '', 'start': 2799, 'end': 4605, 'color': '#9E51F6'}, {'name': 'WPRE', 'type': 'Miscellaneous', 'description': 'Woodchuck hepatitis virus posttranscriptional regulatory element', 'notes': 'Enhances virus stability in packaging cells, leading to higher titer of packaged virus; enhances higher expression of transgenes.', 'start': 4615, 'end': 5204, 'color': '#9B3BDC'}, {'name': "'3' LTR (ΔU3)'", 'type': 'LTR', 'description': "Truncated HIV-1 3' long terminal repeat", 'notes': "Allows packaging of viral RNA into virus; self-inactivates the 5' LTR by a copying mechanism during viral genome integration; contains polyadenylation signal for transcription termination.", 'start': 5275, 'end': 5509, 'color': '#455FCF'}, {'name': 'SV40 poly(A) signal', 'type': 'PolyA_signal', 'description': 'Simian virus 40 early polyadenylation signal', 'notes': 'Allows transcription termination and polyadenylation of mRNA transcribed by Pol II RNA polymerase.', 'start': 5580, 'end': 5702, 'color': '#C37156'}, {'name': 'AMPR', 'type': 'ORF', 'description': 'AmpiciIIin resistance gene', 'notes': 'Allows E. coli to be resistant to ampiciIIin.', 'start': 6669, 'end': 7530, 'color': '#74B4BA'}, {'name': 'ori', 'type': 'Rep_origin', 'description': 'pUC origin of replication', 'notes': 'Facilitates plasmid replication in E. coli; regulates high-copy plasmid number (500-700).', 'start': 7700, 'end': 8289, 'color': '#8515E1'}]
info_dict={'type': '基因编辑质粒载体（单gRNA）', 'marker': 'EGFP，Puro', 'resist': 'Ampicillin', 'host': 'stbl3', 'description': '哺乳动物基因编辑载体，表达gRNA和Cas9元件，另有荧光和抗性基因可用作筛选和标记。'}
# 页眉
def header(canvas, doc):

    canvas.saveState()
    pageNumber = canvas.getPageNumber()
    p = Paragraph("<img src='%s' width='%d' height='%d'/>" % (r'C:\Users\41518\Desktop\work\Ubigene\carrier\title_head.png', 300, 30),style=getSampleStyleSheet()['Normal'])  # 使用一个Paragraph Flowable存放图片
    w, h = p.wrap(doc.width, doc.bottomMargin)
    p.drawOn(canvas, doc.leftMargin, 790)

    canvas.drawString(10 * cm, cm, str(pageNumber))
    canvas.line(doc.leftMargin, doc.bottomMargin+doc.height + 0.5*cm, doc.leftMargin+doc.width, doc.bottomMargin+doc.height + 0.5*cm) #画一条横线
    canvas.restoreState()



def rpt(result_dict,info_dict):
    story=[]
    stylesheet=getSampleStyleSheet()

    # 设置文本格式
    normalStyle = stylesheet['Normal']
    infoStyle = ParagraphStyle(name='info',fontSize=10,fontName='msyh',leading=12,)

    head1 = '''<para fontSize=10><font face="msyh" >载体基本信息</font><br/></para>'''
    head2 = '''<para fontSize=10><font face="msyh" >载体元件</font><br/></para>'''
    head3 = '''<para fontSize=10><font face="msyh" >载体序列</font><br/></para>'''
    # 载体基本信息表格
    info_data = []
    for item in info_dict.items():
        info_data.append([item[0],Paragraph(item[1],infoStyle)])
    info_table = Table(info_data, colWidths=[100, 400])
    info_table.setStyle(TableStyle([
        ('FONTNAME', (0, 0), (-1, -1), 'msyh'),  # 字体
        ('FONTSIZE', (0, 0), (-1, -1), 10),  # 字体大小

        ('ALIGN', (0, 0), (-1, 0), 'LEFT'),  # 对齐
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),  # 对齐
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),  # 设置表格框线为灰，线宽为0.5
    ]))


    # drawing = svg2rlg('C://Users//41518//Desktop//work//Ubigene//carrier//pdf_vector.svg')
    drawing = svg2rlg('C://Users//41518//Desktop//work//Ubigene//carrier//vector.svg')
    drawing.width, drawing.height = drawing.minWidth() * 0.6, drawing.height * 0.6
    drawing.scale(0.6, 0.6)

    # #图片，用法详见reportlab-userguide.pdf中chapter 9.3 Image
    # img = Image('C://Users//41518//Desktop//work//Ubigene//carrier//vector.svg', width=2*inch, height=2*inch)

    # 载体元件表格
    component_data= [['名称', '位置', '大小(bp)', '类型','描述']]
    for item in result_dict:
        content1 = "<para><font color={}>。 {}-{}</font><font></font></para>".format(item['color'],item['start'], item['end'])
        component_data.append([item['name'],Paragraph(content1,normalStyle),item['end']-item['start']+1,item['type'],Paragraph(item['description'],normalStyle)])
    #创建表格对象，并设定各列宽度
    component_table = Table(component_data,colWidths=[100,80,50,110,160])
    #添加表格样式
    component_table.setStyle(TableStyle([
    ('FONTNAME',(0,0),(-1,-1),'msyh'),#字体
    ('FONTSIZE',(0,0),(-1,-1),10),#字体大小
    # ('SPAN',(0,0),(3,0)),#合并第一行前三列
    ('BACKGROUND',(0,0),(-1,0), colors.lightskyblue),#设置第一行背景颜色
    # ('SPAN',(-1,0),(-2,0)), #合并第一行后两列
    ('ALIGN',(0,0),(-1,0),'CENTER'),#对齐
    ('VALIGN',(0,0),(-1,-1),'MIDDLE'),  #对齐
    ('LINEBEFORE',(0,0),(0,-1),0.1,colors.grey),#设置表格左边线颜色为灰色，线宽为0.1
    # ('TEXTCOLOR',(0,1),(-2,-1),colors.royalblue),#设置表格内文字颜色
    ('GRID',(0,0),(-1,-1),0.5,colors.red),#设置表格框线为红色，线宽为0.5
    ]))



    story.append(Paragraph(head1, normalStyle))
    story.append(info_table)
    story.append(drawing)
    story.append(Paragraph(head2,normalStyle))
    story.append(component_table)

    # 设置页边距，页眉
    doc = SimpleDocTemplate('/Ubigene/carrier/双gRNA.pdf', leftMargin=0.6 * cm, rightMargin=0.6 * cm)
    frame_footer = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
    template = PageTemplate(id='story',frames=frame_footer ,onPage=header)
    doc.addPageTemplates([template])
    doc.build(story)

if __name__ == '__main__':
    rpt(result_dict,info_dict)