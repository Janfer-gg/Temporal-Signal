import pymysql
import json
import re

# 完成设计后选定载体
def get_vector(scaffold_type,marker_name,protein_name):
    global result
    if marker_name and protein_name:
        cursor.execute("SELECT vector FROM vectorbuilder1 WHERE type = '%s' and marker = '%s' and protein = '%s'  " % (scaffold_type,marker_name,protein_name))
        result = cursor.fetchone()[0]
    elif marker_name:
        cursor.execute("SELECT vector FROM vectorbuilder1 WHERE type = '%s' and marker = '%s'   " % (scaffold_type, marker_name))
        result = cursor.fetchone()[0]
    elif protein_name:
        cursor.execute("SELECT vector FROM vectorbuilder1 WHERE type = '%s'  and protein = '%s'  " % (scaffold_type, protein_name))
        result = cursor.fetchone()[0]
    return result

# 完成设计后获取载体信息
def get_details(scaffold_type,vector,gRNA_seq1,gRNA_name1,gRNA_seq2='',gRNA_name2=''):
    cursor.execute("SELECT * FROM vectorbuilder1 WHERE type = '%s' and vector = '%s'  " % (scaffold_type,vector))
    result = cursor.fetchone()
    seq_dict = eval(result[5])
    seq = result[6]
    scaffold_info = eval(result[7])
    details = eval(result[8])
    if scaffold_type == '表达慢病毒载体单gRNA' or scaffold_type == 'CRISPR慢病毒载体单gRNA':
        details[6]['name'] = gRNA_name1
        for i,j in zip(seq_dict.items(),details):
            j['start'] = re.search(i[1],seq,flags=re.I).span()[0]+1
            j['end'] = re.search(i[1], seq, flags=re.I).span()[1]

        seq = re.sub('[nN]+', gRNA_seq1, seq)

    elif scaffold_type == '表达慢病毒载体双gRNA' or scaffold_type == 'CRISPR慢病毒载体双gRNA':
        details[6]['name'] = gRNA_name1
        details[8]['name'] = gRNA_name2

        for i, j in zip(seq_dict.items(), details):
            j['start'] = re.search(i[1], seq, flags=re.I).span()[0] + 1
            j['end'] = re.search(i[1], seq, flags=re.I).span()[1]
        seq = re.sub('[nN]+', gRNA_seq1, seq)
        seq = re.sub('[mM]+', gRNA_seq2, seq)

    elif scaffold_type == 'CRISPR载体单gRNA':
        details[1]['name'] = gRNA_name1
        for i,j in zip(seq_dict.items(),details):
            j['start'] = re.search(i[1],seq,flags=re.I).span()[0]+1
            j['end'] = re.search(i[1], seq, flags=re.I).span()[1]

        seq = re.sub('[nN]+', gRNA_seq1, seq)

    elif scaffold_type == 'CRISPR载体双gRNA':
        details[1]['name'] = gRNA_name1
        details[3]['name'] = gRNA_name2
        for i, j in zip(seq_dict.items(), details):
            j['start'] = re.search(i[1], seq, flags=re.I).span()[0] + 1
            j['end'] = re.search(i[1], seq, flags=re.I).span()[1]
        seq = re.sub('[nN]+', gRNA_seq1, seq)
        seq = re.sub('[mM]+', gRNA_seq2, seq)

    return [scaffold_info,details,seq]


if __name__ == '__main__':
    db = pymysql.connect("47.92.215.36", "ubigene_data", "@ubigene2020@..", "YuanJingData", charset='utf8')
    cursor = db.cursor()
    scaffold_type='CRISPR载体双gRNA'
    marker_name='EGFP:P2A:Puro'
    protein_name='hCas9'
    gRNA_name1 = 'sample'
    gRNA_seq1 = 'AAAAAAAAAAAAAAAAAAAAA'
    vector = get_vector(scaffold_type,marker_name,protein_name)
    aa = get_details(scaffold_type,vector,gRNA_seq1,gRNA_name1)
    print(aa)
