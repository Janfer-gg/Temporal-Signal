import re
import math
import random
import pymysql
from reportlab.graphics import renderPDF
from svglib.svglib import svg2rlg

def get_angle(bp, length):
    """质粒碱基的位置转换为角度"""
    return bp * 360 / length

def coord(angle,center, radius):
    #角度转弧度
    rad1 = math.radians(90-angle)
    x = center[0] + math.cos(rad1) * radius
    y = center[1] - math.sin(rad1) * radius
    return x, y

#获取箭头
def arrowpath(start,end,CENTER,radius,PLASMID_LENGTH):
    angle1 = get_angle(start,PLASMID_LENGTH)
    angle2 = get_angle(end,PLASMID_LENGTH)
    angle3 = get_angle(end - 54,PLASMID_LENGTH)
    p1 = coord(angle1, CENTER, radius[1])
    p2 = coord(angle3, CENTER, radius[1])
    p3 = coord(angle3, CENTER, radius[0])
    p4 = coord(angle2, CENTER, radius[3])
    p5 = coord(angle3, CENTER, radius[6])
    p6 = coord(angle3, CENTER, radius[5])
    p7 = coord(angle1, CENTER, radius[5])
    if(end-start>=50):
        d1 = 'M %s %s A %d %d 0 0 1 %s %s L %s %s L %s %s L %s %s L %s %s A %d %d 0 0 0 %s %s L %s %s z' % (p1[0], p1[1], radius[1], radius[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], p5[0], p5[1], p6[0], p6[1],radius[5], radius[5], p7[0], p7[1], p1[0], p1[1])
    else:
        d1='M %s %s L %s %s'%(p1[0], p1[1], p7[0], p7[1])
    return(d1)
#
def text_coord(start1,end1,start2,end2,type,PLASMID_LENGTH):
    if start2 in range(int(PLASMID_LENGTH/8),int(PLASMID_LENGTH*3/8)) or start2 in range(int(PLASMID_LENGTH*5/8),int(PLASMID_LENGTH*7/8)):
        if start1 != 0:
            if type == False:
                if ((start2 + end2) - (start1 + end1) > 200):
                    return False
                else:
                    return True
            if type == True:
                if ((start2 + end2) - (start1 + end1) > 400):
                    return False
                else:
                    return True
        else:
            return False
    else:
        if start1 != 0:
            if type == False:
                if ((start2 + end2) - (start1 + end1) > 500):
                    return False
                else:
                    return True
            if type == True:
                if ((start2 + end2) - (start1 + end1) > 1000):
                    return False
                else:
                    return True
        else:
            return False

#获取文本注释线
def markerline(start, end, type, CENTER, radius, PLASMID_LENGTH):
    angle4 = get_angle((start + end) / 2, PLASMID_LENGTH)
    angle5 = get_angle((start + end) / 2 + 500, PLASMID_LENGTH)
    p1 = coord(angle4, CENTER, radius[5])
    p2 = coord(angle4, CENTER, radius[7])
    if type:
        p3 = coord(angle5, CENTER, radius[8])
    else:
        p3 = coord(angle4, CENTER, radius[8])

    if (start + end) < PLASMID_LENGTH:
        d2 = 'M %s %s L %s %s L %s %s L %s %s' % (p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p3[0] + 5, p3[1])
    else:
        d2 = 'M %s %s L %s %s L %s %s L %s %s' % (p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p3[0] - 5, p3[1])
    return d2

#获取文本坐标
def text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH):
    angle4 = get_angle((start + end) / 2, PLASMID_LENGTH)
    angle5 = get_angle((start + end) / 2 + 500, PLASMID_LENGTH)

    if type:
        p1 = coord(angle5, CENTER, radius[8])
    else:
        p1 = coord(angle4, CENTER, radius[8])

    if (start + end) < PLASMID_LENGTH:
        p2 = (p1[0] + 10, p1[1] + 4)
    else:
        p2 = (p1[0] - 10, p1[1] + 4)
    return p2

#获取文本方向
def text_anchor(start,end,PLASMID_LENGTH):
    if (start+end) < PLASMID_LENGTH:
        text = 'start'
    else:
        text = 'end'
    return text

#获取颜色
def randomcolor(length):
    color_list=[]
    for i in range(length):
        colorArr = ['1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
        color = ""
        for i in range(6):
            color += colorArr[random.randint(0,14)]
        color_list.append("#"+color)
    return color_list

#svg图
def get_svg(scaffold_type,scaffold_name,vector_name):
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8')
    cursor = db.cursor()
    cursor.execute("SELECT * FROM vector_reporter WHERE type = '%s' AND vector = '%s'  " % (scaffold_type, scaffold_name))
    result = cursor.fetchone()
    seq = result[3]
    seq_dict = eval(result[2])
    gRNA_label = re.findall('gRNA\d+|shRNA\d+',vector_name)

    dict_key_list = list(seq_dict.keys())
    dict_value_list = list(seq_dict.values())
    if scaffold_type == 'dual gRNA':
        dict_key_list = [gRNA_label[0] if i == 'gRNA' else i for i in dict_key_list]
        dict_key_list = [gRNA_label[1] if i == 'gRNA ' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))
    else:
        dict_key_list = [gRNA_label[0] if i == 'gRNA' or i == 'shRNA' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))
    PLASMID_LENGTH = seq.__len__()
    CENTER = (400, 250)
    radius = (133, 138, 143, 145.5, 148, 153, 158, 175, 220)
    start = 0
    end = 0
    type = False
    svg_list=[]
    compo_color=randomcolor(len(seq_dict))
    for compo, color in zip(seq_dict.items(), compo_color):
        type = text_coord(start, end, re.search(compo[1], seq,flags=re.I).span()[0], re.search(compo[1], seq,flags=re.I).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq,flags=re.I).span()[0]
        end = re.search(compo[1], seq,flags=re.I).span()[1]

        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "arrow_path_d": compo_arrow, "arrow_path_style": "fill:{};fill-opacity:.5;stroke:{};stroke-with:2;filter:url(#f1)".format(color,color),
                    "markerline_path_d": compo_marker,"markerlabel_text_anchor": compo_text_anchor,
                    "markerlabel_text_start": compo_text[0],"markerlabel_text_end": compo_text[1]}
        # print(svg_dict)
        svg_list.append(svg_dict)
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" class="plasmid ng-isolate-scope" width="800" height="500" ng-class="{'svgBkg': opts.vectorStyle.trackRotate&amp;&amp;rotShow}" >
    <style type="text/css">
           .markerlabel {
               fill: rgba(0, 0, 64, .7);
               cursor: pointer;
               font-size: 12px;
               font-family: Arial;
               stroke-width: 0;
           }

           .markerline {
               stroke: #666;
               stroke-width: 0.7;
               fill: none;
           }

           .scale1 {
               stroke: #999;
           }

           .scale2 {
               font-size: 9px;
               fill: #999;
               font-weight: 300;
               stroke-width: 0
           }

           .plasmid .marker-off {
               stroke: #d61111;
               stroke-width: 0.3px;
           }

           .pla-outRing-s1 {
               fill: #ccc;
               stroke: #999;
           }

           .trackCen {
               font-size: 20px;
               font-weight: 400;
               fill: #000;
               stroke-width: 0;
           }

           .trackCenBp {
               font-size: 10px;
               fill: #000;
               stroke-width: 0;
           }
           tspan.red {
               fill: red;
           }
           tspan.blue {
               fill: blue;
           }
           tspan.bold {
               font-weight: bold;
           }
       </style>
       <defs>
            <filter id="f1" x="-50%" y="-50%" width="200%" height="200%">
                <feOffset result="offOut" in="SourceGraphic"  />
                <feGaussianBlur result="blurOut" in="offOut" stdDeviation="3" />
                <feBlend in="SourceGraphic" in2="blurOut" mode="normal" />
            </filter>
	    </defs>
           <g class="pla-outRing-s1"><path class="ng-scope ng-isolate-scope pla-outRing-s1" fill-rule="evenodd" d="M 399.975 107 A 143 143 0 1 0 400 107 M 399.974 102 A 148 148 0 1 0 400 102" /></g>\n'''

    for item in svg_list:
        svg = svg + "<g><path class='ng-scope ng-isolate-scope' " + '''d="{}"'''.format(
            item['arrow_path_d']) + ''' style="{}"'''.format(
            item['arrow_path_style']) + "></path>\n<g><path" + ''' d="{}" class="markerline"'''.format(
            item['markerline_path_d']) + "></path>" + "<text" + ''' text-anchor="{}" x="{}" y="{}">'''.format(
            item['markerlabel_text_anchor'], item['markerlabel_text_start'], item['markerlabel_text_end']) + item[
                   'name'] + "</text></g></g>\n"
    # 载体名称分两行
    if len(vector_name) > 25:
        next_line = int(re.search('\[.*\]',vector_name).start())
        svg = svg + '''<g><text text-anchor="middle"><tspan x="400" y="250">{}</tspan><tspan x="400" y="270">{}</tspan></text></g>\n</svg>'''.format(vector_name[0:next_line],vector_name[next_line:])
    # 一行
    else:
        svg = svg + '''<g><text text-anchor="middle" x="400" y="250">{}</text></g>\n</svg>'''.format(vector_name)

    with open("pdf_vector.svg", 'w', encoding='utf-8') as f:
        f.write(svg)
        f.close()



if __name__ == '__main__':
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8' )
    cursor = db.cursor()

    scaffold_type = "single gRNA"
    scaffold_name = 'YKO-LV003'
    vector_name = 'YKO-LV003-BCLAF1[gRNA1-gRNA2]'

    get_svg(scaffold_type,scaffold_name,vector_name)




