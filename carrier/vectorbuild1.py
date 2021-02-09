import re
import math
import random
import pymysql
import json


def get_angle(bp, length):
    """质粒碱基的位置转换为角度"""
    return bp * 360 / length


def coord(angle, center, radius):
    # 角度转弧度
    rad1 = math.radians(90 - angle)
    x = center[0] + math.cos(rad1) * radius
    y = center[1] - math.sin(rad1) * radius
    return x, y


# 获取箭头
def arrowpath(start, end, CENTER, radius, PLASMID_LENGTH):
    angle1 = get_angle(start, PLASMID_LENGTH)
    angle2 = get_angle(end, PLASMID_LENGTH)
    angle3 = get_angle(end - 54, PLASMID_LENGTH)
    p1 = coord(angle1, CENTER, radius[1])
    p2 = coord(angle3, CENTER, radius[1])
    p3 = coord(angle3, CENTER, radius[0])
    p4 = coord(angle2, CENTER, radius[2])
    p5 = coord(angle3, CENTER, radius[4])
    p6 = coord(angle3, CENTER, radius[3])
    p7 = coord(angle1, CENTER, radius[3])
    if (end - start >= 60):
        d1 = 'M %s %s A %d %d 0 0 1 %s %s L %s %s L %s %s L %s %s L %s %s A %d %d 0 0 0 %s %s L %s %s z' % (
        p1[0], p1[1], radius[1], radius[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1], p5[0], p5[1], p6[0], p6[1],
        radius[5], radius[5], p7[0], p7[1], p1[0], p1[1])
    else:
        d1 = 'M %s %s L %s %s' % (p1[0], p1[1], p7[0], p7[1])
    return (d1)


#
def text_coord(start1, end1, start2, end2, type, PLASMID_LENGTH):
    if start2 in range(int(PLASMID_LENGTH / 8), int(PLASMID_LENGTH * 3 / 8)) or start2 in range(
            int(PLASMID_LENGTH * 5 / 8), int(PLASMID_LENGTH * 7 / 8)):
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

# 获取文本注释线
def markerline(start, end, type, CENTER, radius, PLASMID_LENGTH):
    angle4 = get_angle((start + end) / 2, PLASMID_LENGTH)
    angle5 = get_angle((start + end) / 2 + 200, PLASMID_LENGTH)
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
# 获取文本注释线
def markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH):
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

# 获取文本坐标
def text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH):
    angle4 = get_angle((start + end) / 2, PLASMID_LENGTH)
    angle5 = get_angle((start + end) / 2 + 200, PLASMID_LENGTH)

    if type:
        p1 = coord(angle5, CENTER, radius[8])
    else:
        p1 = coord(angle4, CENTER, radius[8])

    if (start + end) < PLASMID_LENGTH:
        p2 = (p1[0] + 10, p1[1] + 4)
    else:
        p2 = (p1[0] - 10, p1[1] + 4)
    return p2
# 获取文本坐标
def text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH):
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


# 获取文本方向
def text_anchor(start, end, PLASMID_LENGTH):
    if (start + end) < PLASMID_LENGTH:
        text = 'start'
    else:
        text = 'end'
    return text


# 获取颜色
def randomcolor(length):
    color_list = []
    for i in range(length):
        colorArr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
        color = ""
        for i in range(6):
            color += colorArr[random.randint(0, 14)]
        color_list.append("#" + color)
    return color_list

# 获取元件序列
def get_component(type, name):
    if type:
        cursor.execute("SELECT `详细信息` FROM vectorbuilder WHERE `种类` = '%s'  AND `元件` = '%s'  " % (type, name))
    else:
        cursor.execute("SELECT `详细信息` FROM vectorbuilder WHERE `元件` = '%s'  " % (name))
    result = re.sub(rf"(true|false|null)", r" '\1'", cursor.fetchone()[0])
    result_dict = eval(result)
    component_seq = result_dict['seq']
    return component_seq


# 表达慢病毒载体单gRNA
def get_svg(scaffold_type, gRNA_name, gRNA_seq, marker_type, marker_name):
    cursor.execute("SELECT * FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    result = cursor.fetchone()
    seq = result[1]
    seq_dict = eval(result[2])

    if gRNA_name and gRNA_seq:
        seq = re.sub('[N]+', gRNA_seq, seq)
        seq_dict['添加 gRNA'] = re.sub('[N]+', gRNA_seq, seq_dict['添加 gRNA'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name if i == '添加 gRNA' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if marker_name and marker_type:
        marker_seq = get_component(marker_type, marker_name)
        seq = re.sub('M', marker_seq, seq)
        seq_dict['添加元件 marker'] = marker_seq
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [marker_name if i == '添加元件 marker' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    PLASMID_LENGTH = seq.__len__()
    CENTER = (500, 300)
    radius = (124,128,131.5,135, 139, 135, 148, 163, 192)
    start = 0
    end = 0
    type = False
    svg_list = [{"sequencelength": PLASMID_LENGTH}]
    compo_color = randomcolor(len(seq_dict))

    for compo, color in zip(seq_dict.items(), compo_color):
        type = text_coord(start, end, re.search(compo[1], seq).span()[0], re.search(compo[1], seq).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq).span()[0]
        end = re.search(compo[1], seq).span()[1]

        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "opacity: 0.5;fill:" + color,
                    "arrow_path_d": compo_arrow, "arrow_path_style": "opacity: 0.8;fill:" + color,
                    "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                    "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                    "markerlabel_text_end": compo_text[1]}
        # print(svg_dict)

        svg_list.append(svg_dict)

    if gRNA_name:
        svg_list[7]['markerlabel_text_style'] = 'fill:' + compo_color[6] + ';font-weight:800'
    else:
        svg_list[7]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name:
        svg_list[9]['markerlabel_text_style'] = 'fill:' + compo_color[8] + ';font-weight:800'
    else:
        svg_list[9]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name and gRNA_name:
        cursor.execute(
            "SELECT * FROM vectorbuilder1 WHERE type = '%s'  AND marker = '%s' " % (scaffold_type, marker_name))
        result = cursor.fetchone()
        seq = re.sub('[n]+', gRNA_seq, result[6])
        PLASMID_LENGTH = seq.__len__()
        CENTER = (500, 300)
        radius = (113, 118, 123, 125.5, 128, 133, 138, 143, 172)
        start = 0
        end = 0
        type = False
        svg_list2 = [{"sequencelength": PLASMID_LENGTH}]
        for compo, color in zip(seq_dict.items(), compo_color):
            type = text_coord(start, end, re.search(compo[1], seq, flags=re.I).span()[0],
                              re.search(compo[1], seq, flags=re.I).span()[1], type,
                              PLASMID_LENGTH)

            start = re.search(compo[1], seq, flags=re.I).span()[0]
            end = re.search(compo[1], seq, flags=re.I).span()[1]

            compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
            compo_marker = markerline(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text = text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

            svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "fill:" + color,
                        "arrow_path_d": compo_arrow, "arrow_path_style": "fill:" + color,
                        "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                        "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                        "markerlabel_text_end": compo_text[1]}

            svg_list2.append(svg_dict)

        return [svg_list, svg_list2]
    else:
        return [svg_list]

# 表达慢病毒载体双gRNA
def get_svg2(scaffold_type, gRNA_name1, gRNA_name2, gRNA_seq1, gRNA_seq2, marker_type, marker_name):
    cursor.execute("SELECT * FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    result = cursor.fetchone()
    seq = result[1]
    seq_dict = eval(result[2])

    if gRNA_name1 and gRNA_seq1:
        seq = re.sub('[N]+', gRNA_seq1, seq)
        seq_dict['添加 gRNA1'] = re.sub('[N]+', gRNA_seq1, seq_dict['添加 gRNA1'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name1 if i == '添加 gRNA1' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if gRNA_name2 and gRNA_seq2:
        seq = re.sub('[H]+', gRNA_seq2, seq)
        seq_dict['添加 gRNA2'] = re.sub('[H]+', gRNA_seq2, seq_dict['添加 gRNA2'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name2 if i == '添加 gRNA2' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if marker_name and marker_type:
        marker_seq = get_component(marker_type, marker_name)
        seq = re.sub('M', marker_seq, seq)
        seq_dict['添加元件 marker'] = marker_seq
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [marker_name if i == '添加元件 marker' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    PLASMID_LENGTH = seq.__len__()
    CENTER = (500, 300)
    radius = (133, 138, 143, 145.5, 148, 153, 158, 163, 192)
    start = 0
    end = 0
    type = False
    svg_list = [{"sequencelength": PLASMID_LENGTH}]
    compo_color = randomcolor(len(seq_dict))
    print(seq_dict)
    for compo, color in zip(seq_dict.items(), compo_color):

        type = text_coord(start, end, re.search(compo[1], seq).span()[0], re.search(compo[1], seq).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq).span()[0]
        end = re.search(compo[1], seq).span()[1]
        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "opacity: 0.5;fill:" + color,
                    "arrow_path_d": compo_arrow, "arrow_path_style": "opacity: 0.5;fill:" + color,
                    "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                    "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                    "markerlabel_text_end": compo_text[1]}
        if compo[0] == "U6 promoter2":
            svg_dict['name'] = "U6 promoter"
        print(svg_dict)
        svg_list.append(svg_dict)

    if gRNA_name1:
        svg_list[7]['markerlabel_text_style'] = 'fill:' + compo_color[6] + ';font-weight:800'
    else:
        svg_list[7]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if gRNA_name2:
        svg_list[9]['markerlabel_text_style'] = 'fill:' + compo_color[8] + ';font-weight:800'
    else:
        svg_list[9]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name:
        svg_list[11]['markerlabel_text_style'] = 'fill:' + compo_color[10] + ';font-weight:800'
    else:
        svg_list[11]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[11]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[11]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name and gRNA_name1 and gRNA_name2:
        cursor.execute("SELECT * FROM vectorbuilder1 WHERE type = '%s'  AND marker = '%s' " % (scaffold_type, marker_name))
        result = cursor.fetchone()
        seq = re.sub('[n]+', gRNA_seq1, result[6])
        seq = re.sub('[m]+', gRNA_seq2, seq)
        PLASMID_LENGTH = seq.__len__()
        CENTER = (500, 300)
        radius = (113, 118, 123, 125.5, 128, 133, 138, 143, 172)
        start = 0
        end = 0
        type = False
        svg_list2 = [{"sequencelength": PLASMID_LENGTH}]
        for compo, color in zip(seq_dict.items(), compo_color):
            type = text_coord(start, end, re.search(compo[1], seq, flags=re.I).span()[0],
                              re.search(compo[1], seq, flags=re.I).span()[1], type,
                              PLASMID_LENGTH)

            start = re.search(compo[1], seq, flags=re.I).span()[0]
            end = re.search(compo[1], seq, flags=re.I).span()[1]

            compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
            compo_marker = markerline(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text = text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

            svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "fill:" + color,
                        "arrow_path_d": compo_arrow, "arrow_path_style": "fill:" + color,
                        "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                        "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                        "markerlabel_text_end": compo_text[1]}

            if compo[0] == "U6 promoter2":
                svg_dict["name"] = "U6 promoter"

            svg_list2.append(svg_dict)

        return [svg_list, svg_list2]
    else:

        return [svg_list]

# CRISPR慢病毒载体单gRNA
def get_svg3(scaffold_type, gRNA_name, gRNA_seq, protein_type, protein_name):
    cursor.execute("SELECT * FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    result = cursor.fetchone()
    seq = result[1]
    seq_dict = eval(result[2])

    if gRNA_name and gRNA_seq:
        seq = re.sub('[N]+', gRNA_seq, seq)
        seq_dict['添加 gRNA'] = re.sub('[N]+', gRNA_seq, seq_dict['添加 gRNA'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name if i == '添加 gRNA' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if protein_name and protein_type:
        protein_seq = get_component(protein_type, protein_name)
        seq = re.sub('[M]+', protein_seq, seq)
        seq_dict['添加元件 Cas protein'] = protein_seq
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [protein_name if i == '添加元件 Cas protein' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    PLASMID_LENGTH = seq.__len__()
    CENTER = (500, 300)
    radius = (124,128,131.5,135, 139, 143, 148, 163, 192)
    start = 0
    end = 0
    type = False
    svg_list = [{"sequencelength": PLASMID_LENGTH}]
    compo_color = randomcolor(len(seq_dict))

    for compo, color in zip(seq_dict.items(), compo_color):
        type = text_coord(start, end, re.search(compo[1], seq).span()[0], re.search(compo[1], seq).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq).span()[0]
        end = re.search(compo[1], seq).span()[1]

        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "opacity: 0.5;fill:" + color,
                    "arrow_path_d": compo_arrow, "arrow_path_style": "opacity: 0.5;fill:" + color,
                    "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                    "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                    "markerlabel_text_end": compo_text[1]}
        # print(svg_dict)
        svg_list.append(svg_dict)

    if gRNA_name:
        svg_list[7]['markerlabel_text_style'] = 'fill:' + compo_color[6] + ';font-weight:800'
    else:
        svg_list[7]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if protein_name:
        svg_list[9]['markerlabel_text_style'] = 'fill:' + compo_color[8] + ';font-weight:800'
    else:
        svg_list[9]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    # 小图
    if protein_name and gRNA_name:
        cursor.execute("SELECT * FROM vectorbuilder1 WHERE type = '%s'  AND protein = '%s' " % (scaffold_type, protein_name))
        result = cursor.fetchone()
        seq = re.sub('[n]+', gRNA_seq, result[6])
        PLASMID_LENGTH = seq.__len__()
        CENTER = (500, 300)
        radius = (113, 118, 123, 125.5, 128, 133, 138, 143, 200)
        start = 0
        end = 0
        type = False
        svg_list2 = [{"sequencelength": PLASMID_LENGTH}]

        for compo, color in zip(seq_dict.items(), compo_color):
            type = text_coord(start, end, re.search(compo[1], seq, flags=re.I).span()[0],
                              re.search(compo[1], seq, flags=re.I).span()[1], type,
                              PLASMID_LENGTH)

            start = re.search(compo[1], seq, flags=re.I).span()[0]
            end = re.search(compo[1], seq, flags=re.I).span()[1]

            compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
            compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

            svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "fill:" + color,
                        "arrow_path_d": compo_arrow, "arrow_path_style": "fill:" + color,
                        "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                        "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                        "markerlabel_text_end": compo_text[1]}
            # print(svg_dict)
            svg_list2.append(svg_dict)
        return [svg_list, svg_list2]
    else:
        return [svg_list]

# CRISPR慢病毒载体双gRNA
def get_svg4(scaffold_type, gRNA_name1, gRNA_name2, gRNA_seq1, gRNA_seq2, protein_type, protein_name):
    cursor.execute("SELECT * FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    result = cursor.fetchone()
    seq = result[1]
    seq_dict = eval(result[2])

    if gRNA_name1 and gRNA_seq1:
        seq = re.sub('[N]+', gRNA_seq1, seq)
        seq_dict['添加 gRNA1'] = re.sub('[N]+', gRNA_seq1, seq_dict['添加 gRNA1'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name1 if i == '添加 gRNA1' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if gRNA_name2 and gRNA_seq2:
        seq = re.sub('[H]+', gRNA_seq2, seq)
        seq_dict['添加 gRNA2'] = re.sub('[H]+', gRNA_seq2, seq_dict['添加 gRNA2'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name2 if i == '添加 gRNA2' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if protein_name and protein_type:
        protein_seq = get_component(protein_type, protein_name)
        seq = re.sub('[M]+', protein_seq, seq)
        seq_dict['添加元件 Cas protein'] = protein_seq
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [protein_name if i == '添加元件 Cas protein' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    PLASMID_LENGTH = seq.__len__()
    CENTER = (500, 300)
    radius = (124,128,131.5,135, 139, 143, 148, 163, 192)
    start = 0
    end = 0
    type = False
    svg_list = [{"sequencelength": PLASMID_LENGTH}]
    compo_color = randomcolor(len(seq_dict))

    for compo, color in zip(seq_dict.items(), compo_color):

        type = text_coord(start, end, re.search(compo[1], seq).span()[0], re.search(compo[1], seq).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq).span()[0]
        end = re.search(compo[1], seq).span()[1]
        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "opacity: 0.5;fill:" + color,
                    "arrow_path_d": compo_arrow, "arrow_path_style": "opacity: 0.5;fill:" + color,
                    "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                    "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                    "markerlabel_text_end": compo_text[1]}
        if compo[0] == "U6 promoter2":
            svg_dict['name'] = "U6 promoter"
        svg_list.append(svg_dict)

    if gRNA_name1:
        svg_list[7]['markerlabel_text_style'] = 'fill:' + compo_color[6] + ';font-weight:800'
    else:
        svg_list[7]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[7]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if gRNA_name2:
        svg_list[9]['markerlabel_text_style'] = 'fill:' + compo_color[8] + ';font-weight:800'
    else:
        svg_list[9]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[9]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if protein_name:
        svg_list[11]['markerlabel_text_style'] = 'fill:' + compo_color[10] + ';font-weight:800'
    else:
        svg_list[11]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[11]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[11]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if protein_name and gRNA_name1 and gRNA_name2:
        cursor.execute(
            "SELECT * FROM vectorbuilder1 WHERE type = '%s'  AND protein = '%s' " % (scaffold_type, protein_name))
        result = cursor.fetchone()
        seq = re.sub('[n]+', gRNA_seq1, result[6])
        seq = re.sub('[m]+', gRNA_seq2, seq)
        PLASMID_LENGTH = seq.__len__()
        CENTER = (500, 300)
        radius = (113, 118, 123, 125.5, 128, 133, 138, 143, 200)
        start = 0
        end = 0
        type = False
        svg_list2 = [{"sequencelength": PLASMID_LENGTH}]

        for compo, color in zip(seq_dict.items(), compo_color):
            type = text_coord(start, end, re.search(compo[1], seq, flags=re.I).span()[0],
                              re.search(compo[1], seq, flags=re.I).span()[1], type,
                              PLASMID_LENGTH)

            start = re.search(compo[1], seq, flags=re.I).span()[0]
            end = re.search(compo[1], seq, flags=re.I).span()[1]

            compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
            compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

            svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "fill:" + color,
                        "arrow_path_d": compo_arrow, "arrow_path_style": "fill:" + color,
                        "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                        "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                        "markerlabel_text_end": compo_text[1]}

            if compo[0] == "U6 promoter2":
                svg_dict["name"] = "U6 promoter"

            svg_list2.append(svg_dict)

        return [svg_list, svg_list2]
    else:

        return [svg_list]

# CRISPR载体单gRNA
def get_svg5(scaffold_type, gRNA_name, gRNA_seq, protein_type,protein_name,marker_type, marker_name):
    cursor.execute("SELECT * FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    result = cursor.fetchone()
    seq = result[1]
    seq_dict = eval(result[2])

    if gRNA_name and gRNA_seq:
        seq = re.sub('[N]+', gRNA_seq, seq)
        seq_dict['添加 gRNA'] = re.sub('[N]+', gRNA_seq, seq_dict['添加 gRNA'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name if i == '添加 gRNA' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if protein_type and protein_name:
        protein_seq = get_component(protein_type, protein_name)
        seq = re.sub('[M]+', protein_seq, seq)
        seq_dict['添加元件 Cas protein'] = protein_seq
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [protein_name if i == '添加元件 Cas protein' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))


    if marker_type and marker_name:
        if marker_type == 'No marker' and marker_name == 'No marker':
            dict_key_list = list(seq_dict.keys())
            dict_value_list = list(seq_dict.values())
            dict_key_list = [marker_name if i == '添加元件 marker' else i for i in dict_key_list]
        else:
            marker_seq = get_component(marker_type, marker_name)
            CMV_seq = get_component('','CMV')
            SV40_seq = 'taagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagtt'
            seq = re.sub('K',CMV_seq+'cggatccgccacc'+marker_seq+'ggtacccagacatga'+SV40_seq, seq)
            seq_dict['添加元件 marker'] = marker_seq
            dict_key_list = list(seq_dict.keys())
            dict_value_list = list(seq_dict.values())
            dict_key_list = [marker_name if i == '添加元件 marker' else i for i in dict_key_list]
            dict_key_list.insert(6,'SV40 poly(A) signal')
            dict_value_list.insert(6, SV40_seq)
            dict_key_list.insert(5,'CMV promoter')
            dict_value_list.insert(5, CMV_seq)

        seq_dict = dict(zip(dict_key_list, dict_value_list))


    PLASMID_LENGTH = seq.__len__()
    CENTER = (500, 300)
    radius = (124,128,131.5,135, 139, 143, 148, 163, 192)
    start = 0
    end = 0
    type = False
    svg_list = [{"sequencelength": PLASMID_LENGTH}]
    compo_color = randomcolor(len(seq_dict))

    for compo, color in zip(seq_dict.items(), compo_color):
        type = text_coord(start, end, re.search(compo[1], seq).span()[0], re.search(compo[1], seq).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq).span()[0]
        end = re.search(compo[1], seq).span()[1]

        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "opacity: 0.5;fill:" + color,
                    "arrow_path_d": compo_arrow, "arrow_path_style": "opacity: 0.5;fill:" + color,
                    "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                    "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                    "markerlabel_text_end": compo_text[1]}
        # print(svg_dict)
        svg_list.append(svg_dict)

    if gRNA_name:
        svg_list[2]['markerlabel_text_style'] = 'fill:' + compo_color[1] + ';font-weight:800'
    else:
        svg_list[2]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[2]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[2]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if protein_name:
        svg_list[4]['markerlabel_text_style'] = 'fill:' + compo_color[3] + ';font-weight:800'
    else:
        svg_list[4]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[4]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[4]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name:
        svg_list[7]['markerlabel_text_style'] = 'fill:' + compo_color[6] + ';font-weight:800'
    else:
        svg_list[6]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[6]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[6]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    # 完成设计
    if protein_name and gRNA_name and marker_name:
        cursor.execute("SELECT * FROM vectorbuilder1 WHERE type = '%s'  AND marker = '%s' AND protein = '%s' " % (scaffold_type, marker_name,protein_name))
        result = cursor.fetchone()
        seq = re.sub('[n]+', gRNA_seq, result[6])
        PLASMID_LENGTH = seq.__len__()
        CENTER = (500, 300)
        radius = (113, 118, 123, 125.5, 128, 133, 138, 143, 172)
        start = 0
        end = 0
        type = False
        svg_list2 = [{"sequencelength": PLASMID_LENGTH}]
        if marker_name == 'No marker':
            del seq_dict['No marker']
        for compo, color in zip(seq_dict.items(), compo_color):
            type = text_coord(start, end, re.search(compo[1], seq, flags=re.I).span()[0],
                              re.search(compo[1], seq, flags=re.I).span()[1], type,
                              PLASMID_LENGTH)

            start = re.search(compo[1], seq, flags=re.I).span()[0]
            end = re.search(compo[1], seq, flags=re.I).span()[1]

            compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
            compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

            svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "fill:" + color,
                        "arrow_path_d": compo_arrow, "arrow_path_style": "fill:" + color,
                        "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                        "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                        "markerlabel_text_end": compo_text[1]}
            # print(svg_dict)
            svg_list2.append(svg_dict)
        return [svg_list, svg_list2]
    else:
        return [svg_list]

# CRISPR载体双gRNA
def get_svg6(scaffold_type, gRNA_name1, gRNA_name2, gRNA_seq1, gRNA_seq2, protein_type, protein_name, marker_type, marker_name):
    cursor.execute("SELECT * FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    result = cursor.fetchone()
    seq = result[1]
    seq_dict = eval(result[2])

    if gRNA_name1 and gRNA_seq1:
        seq = re.sub('[N]+', gRNA_seq1, seq)
        seq_dict['添加 gRNA1'] = re.sub('[N]+', gRNA_seq1, seq_dict['添加 gRNA1'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name1 if i == '添加 gRNA1' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if gRNA_name2 and gRNA_seq2:
        seq = re.sub('[H]+', gRNA_seq2, seq)
        seq_dict['添加 gRNA2'] = re.sub('[H]+', gRNA_seq2, seq_dict['添加 gRNA2'])
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [gRNA_name2 if i == '添加 gRNA2' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if protein_type and protein_name:
        protein_seq = get_component(protein_type, protein_name)
        seq = re.sub('[M]+', protein_seq, seq)
        seq_dict['添加元件 Cas protein'] = protein_seq
        dict_key_list = list(seq_dict.keys())
        dict_value_list = list(seq_dict.values())
        dict_key_list = [protein_name if i == '添加元件 Cas protein' else i for i in dict_key_list]
        seq_dict = dict(zip(dict_key_list, dict_value_list))

    if marker_type and marker_name:
        if marker_type == 'No marker' and marker_name == 'No marker':
            dict_key_list = list(seq_dict.keys())
            dict_value_list = list(seq_dict.values())
            dict_key_list = [marker_name if i == '添加元件 marker' else i for i in dict_key_list]
        else:
            marker_seq = get_component(marker_type, marker_name)
            CMV_seq = get_component('', 'CMV')
            SV40_seq = 'taagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagtt'
            seq = re.sub('K', CMV_seq + 'cggatccgccacc' + marker_seq + 'ggtacccagacatga' + SV40_seq, seq)
            seq_dict['添加元件 marker'] = marker_seq
            dict_key_list = list(seq_dict.keys())
            dict_value_list = list(seq_dict.values())
            dict_key_list = [marker_name if i == '添加元件 marker' else i for i in dict_key_list]
            dict_key_list.insert(8, 'SV40 poly(A) signal')
            dict_value_list.insert(8, SV40_seq)
            dict_key_list.insert(7, 'CMV promoter')
            dict_value_list.insert(7, CMV_seq)

        seq_dict = dict(zip(dict_key_list, dict_value_list))
    PLASMID_LENGTH = seq.__len__()
    CENTER = (500, 300)
    radius = (124,128,131.5,135, 139, 143, 148, 163, 192)
    start = 0
    end = 0
    type = False
    svg_list = [{"sequencelength": PLASMID_LENGTH}]
    compo_color = randomcolor(len(seq_dict))

    for compo, color in zip(seq_dict.items(), compo_color):

        type = text_coord(start, end, re.search(compo[1], seq,flags=re.I).span()[0], re.search(compo[1], seq,flags=re.I).span()[1], type,
                          PLASMID_LENGTH)

        start = re.search(compo[1], seq,flags=re.I).span()[0]
        end = re.search(compo[1], seq,flags=re.I).span()[1]
        compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
        compo_marker = markerline2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text = text_pos2(start, end, type, CENTER, radius, PLASMID_LENGTH)
        compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

        svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "opacity: 0.5;fill:" + color,
                    "arrow_path_d": compo_arrow, "arrow_path_style": "opacity: 0.5;fill:" + color,
                    "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                    "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                    "markerlabel_text_end": compo_text[1]}
        if compo[0] == "U6 promoter2":
            svg_dict['name'] = "U6 promoter"
        svg_list.append(svg_dict)

    if gRNA_name1:
        svg_list[2]['markerlabel_text_style'] = 'fill:' + compo_color[1] + ';font-weight:800'
    else:
        svg_list[2]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[2]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[2]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if gRNA_name2:
        svg_list[4]['markerlabel_text_style'] = 'fill:' + compo_color[3] + ';font-weight:800'
    else:
        svg_list[4]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[4]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[4]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if protein_name:
        svg_list[6]['markerlabel_text_style'] = 'fill:' + compo_color[5] + ';font-weight:800'
    else:
        svg_list[6]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[6]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[6]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name:
        svg_list[9]['markerlabel_text_style'] = 'fill:' + compo_color[8] + ';font-weight:800'
    else:
        svg_list[8]['markerstyle'] = 'fill:#fff;stroke:#e04993'
        svg_list[8]['arrow_path_style'] = 'fill:#fff;stroke:#e04993'
        svg_list[8]['markerlabel_text_style'] = 'fill:#e04993;font-weight:700;'

    if marker_name and gRNA_name1 and gRNA_name2 and protein_name:
        cursor.execute("SELECT * FROM vectorbuilder1 WHERE type = '%s'  AND marker = '%s' AND protein = '%s' " % (scaffold_type, marker_name,protein_name))
        result = cursor.fetchone()
        seq = re.sub('[n]+', gRNA_seq1, result[6])
        seq = re.sub('[m]+', gRNA_seq2, seq)
        PLASMID_LENGTH = seq.__len__()
        CENTER = (500, 300)
        radius = (113, 118, 123, 125.5, 128, 133, 138, 143, 200)
        start = 0
        end = 0
        type = False
        svg_list2 = [{"sequencelength": PLASMID_LENGTH}]
        if marker_name == 'No marker':
            del seq_dict['No marker']
        for compo, color in zip(seq_dict.items(), compo_color):
            type = text_coord(start, end, re.search(compo[1], seq, flags=re.I).span()[0],
                              re.search(compo[1], seq, flags=re.I).span()[1], type,
                              PLASMID_LENGTH)

            start = re.search(compo[1], seq, flags=re.I).span()[0]
            end = re.search(compo[1], seq, flags=re.I).span()[1]

            compo_arrow = arrowpath(start, end, CENTER, radius, PLASMID_LENGTH)
            compo_marker = markerline(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text = text_pos(start, end, type, CENTER, radius, PLASMID_LENGTH)
            compo_text_anchor = text_anchor(start, end, PLASMID_LENGTH)

            svg_dict = {"name": compo[0], "start": start, "end": end, "markerstyle": "fill:" + color,
                        "arrow_path_d": compo_arrow, "arrow_path_style": "fill:" + color,
                        "markerlabel_text": compo[0], "markerline_path_d": compo_marker,
                        "markerlabel_text_anchor": compo_text_anchor, "markerlabel_text_start": compo_text[0],
                        "markerlabel_text_end": compo_text[1]}

            if compo[0] == "U6 promoter2":
                svg_dict["name"] = "U6 promoter"

            svg_list2.append(svg_dict)

        return [svg_list, svg_list2]
    else:

        return [svg_list]

if __name__ == '__main__':
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8')
    cursor = db.cursor()

    scaffold_type = "CRISPR载体双gRNA"
    gRNA_name1 = "gRNA"
    gRNA_seq1 = "AAAAAAAAAAAAAAAAAAAA"
    gRNA_name2 = "Sample[gRNA2]"
    gRNA_seq2 = "TTTTTTTTTTTTTTTTTTTT"
    # protein_type = ''
    # protein_name = ''
    # marker_type = ''
    # marker_name = ''
    protein_type = 'Cas Proteins'
    protein_name = 'hCas9'
    marker_type = 'Dual Reporters'
    marker_name = 'EGFP:P2A:Puro'


    if scaffold_type == '表达慢病毒载体单gRNA':
        svg_list = get_svg(scaffold_type, gRNA_name1, gRNA_seq1, marker_type, marker_name)
    elif scaffold_type == '表达慢病毒载体双gRNA':
        svg_list = get_svg2(scaffold_type, gRNA_name1, gRNA_name2, gRNA_seq1, gRNA_seq2, marker_type, marker_name)
    elif scaffold_type == 'CRISPR慢病毒载体单gRNA':
        svg_list = get_svg3(scaffold_type, gRNA_name1, gRNA_seq1, protein_type, protein_name)
    elif scaffold_type == 'CRISPR慢病毒载体双gRNA':
        svg_list = get_svg4(scaffold_type, gRNA_name1, gRNA_name2, gRNA_seq1, gRNA_seq2, protein_type, protein_name)
    elif scaffold_type == 'CRISPR载体单gRNA':
        svg_list = get_svg5(scaffold_type, gRNA_name1,  gRNA_seq1, protein_type, protein_name, marker_type, marker_name)
    elif scaffold_type == 'CRISPR载体双gRNA':
        svg_list = get_svg6(scaffold_type, gRNA_name1, gRNA_name2, gRNA_seq1, gRNA_seq2, protein_type,protein_name, marker_type, marker_name)
    # SVG图形信息
    svg_json = json.dumps(svg_list)
    # print(svg_json)
    # 完成设计后，调取元件信息
    # if gRNA_name1 and marker_name:
    #     cursor.execute("SELECT details FROM vectorbuilder2 WHERE type = '%s'   " % (scaffold_type))
    #     result = cursor.fetchone()
    #     result = result[0]
    #     result_dict = eval(result)
    #
    #     result_dict[6]['name'] = gRNA_name1
    #     result_dict[8]['name'] = marker_name
    #
    #     # print(result_dict)
    #     for i in (range(1,len(svg_list[1]))):
    #         result_dict[i - 1]['start']= svg_list[1][i]['start']
    #         result_dict[i - 1]['end'] = svg_list[1][i]['end']
    #         result_dict[i - 1]['color'] = re.search("\#\w+",svg_list[1][i]['markerstyle']).group()
    #
    # # 元件详细信息
    # print(result_dict)

    list = '''<svg xmlns="http://www.w3.org/2000/svg" class="plasmid ng-isolate-scope" width="1000" height="600" ng-class="{'svgBkg': opts.vectorStyle.trackRotate&amp;&amp;rotShow}" >
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
            <g class="pla-outRing-s1">
                <path class="ng-scope ng-isolate-scope pla-outRing-s1" fill-rule="evenodd" d="M 499.975 157 A 143 143 0 1 0 500 157 M 499.974 152 A 148 148 0 1 0 500 152" />

            </g>\n'''

    for item in svg_list[0]:
        if "sequencelength" in item:
            continue

        list = list + "<g><path class='marker-off markerhover ng-scope ng-isolate-scope' " + '''d="{}"'''.format(item['arrow_path_d']) + ''' style="{}"'''.format(
            item['arrow_path_style']) + "></path>\n<g><path" + ''' d="{}" class="markerline"'''.format(
            item['markerline_path_d']) + "></path>" + "<text " + '''class="markerlabel" text-anchor="{}" x="{}" y="{}">'''.format(
            item['markerlabel_text_anchor'], item['markerlabel_text_start'], item['markerlabel_text_end']) + item[
                   'name'] + "</text></g></g>\n"
    list = list + "</svg>"

    with open("vector.svg", 'w', encoding='utf-8') as f:
        f.write(list)
        f.close()

