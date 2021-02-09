from PIL import Image, ImageFont, ImageDraw
import re

def get_label(type,vector_name,content,size,storage):

    im = Image.open("label.png")  # 打开文件
    im = im.convert('RGB')

    draw = ImageDraw.Draw(im) #修改图片
    font = ImageFont.truetype(r'C:\Users\41518\Desktop\work\Ubigene\标签\msyhBd.ttc', size = 48)
    fillcolor = "black"
    if(len(vector_name)>28):
        gRNA_pos = re.search('\[[\s\S]+\]',vector_name).span()
        vector_name=vector_name[0:gRNA_pos[0]]+'\n'+vector_name[gRNA_pos[0]:gRNA_pos[1]]
        # draw.text((66, 195), vector_name, font=font, fill=fillcolor)
        draw.text((66, 195), vector_name, font=font, fill=fillcolor)
    else:
        draw.text((66,245),vector_name , font=font,fill=fillcolor)

    draw.text((66, 350), 'Content：'+content, font=font, fill=fillcolor)

    if type == 'Plasmid':
        draw.text((66, 430), 'Concentration：'+size, font=font, fill=fillcolor)
    elif type == 'Bacteria':
        draw.text((66, 430), 'Size：' + size, font=font, fill=fillcolor)

    draw.text((66, 510), 'Storage：'+storage, font=font, fill=fillcolor)

    im.save('label2.png')

if __name__ == '__main__':
    vector_name="YKO-RP003[gRNA1-gRNA2]"
    content = 'E.coli，Amp+'
    size = '1ml/Tube'
    storage = '-80℃'

    type = 'Bacteria'
    # content = 'plasmid，Amp+'
    # size = '500ng/ul'
    # storage = '-20℃'

    get_label(type,vector_name,content,size,storage)