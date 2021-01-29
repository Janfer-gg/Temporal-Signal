import os
from copy import copy
from Bio import SeqIO,Align
from Ubigene.ice_analysis.sequence import check_nucleotide,DNA_rev_complement
import warnings
# 输入处理
def input_check(control_path, sample_path, base_outputname, gRNA):
    # 文件处理
    if control_path is None or not os.path.exists(control_path):
        return False,'{} not found'.format(control_path)
        # raise Exception('Control @ {} not found'.format(control_path))

    if os.path.splitext(control_path)[1] not in ['.abi','.ab1']:
        print(os.path.splitext(control_path)[1])
        return False, 'control_file not abi file'
        # raise Exception('control_file not abi file')

    if sample_path is None or not os.path.exists(sample_path):
        return False, '{} not found'.format(sample_path)
        # raise Exception('Experiment sample @ {} not found'.format(sample_path))

    if os.path.splitext(control_path)[1] not in ['.abi','.ab1']:
        return False, 'experiment_file not abi file'
        # raise Exception('sample_file not abi file')

    if base_outputname:
        output_dir = os.path.abspath(base_outputname)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        return False,'Input your output path'
        # raise Exception('Input your output path')

    if gRNA:
        for g in gRNA:
            if 17<=len(g)<=23 and check_nucleotide(g):
                pass
            else:
                return False,'Please check your input gRNA'

        # raise Exception('Input your guide RNA')
    else:
        return False, 'Input your guide RNA'

    return True,''

# abi文件数据获取
def get_data(ab1_file):
    record = SeqIO.read(ab1_file, 'abi')
    data = record.annotations['abif_raw']

    # 质量得分
    quality_scores =[]
    for quality in data['PCON2']:
        quality_scores.append(quality)
    # 碱基顺序
    base_order = str(data['FWO_1'],'utf-8')
    # 主峰序列
    primary_seq = str(data['PBAS1'],'utf-8')
    # 峰的位置
    peak_location =data['PLOC1']
    # 峰值
    output_dict = {}
    for base,channel in zip(base_order,['DATA9','DATA10','DATA11','DATA12']):
        peak_data = data[channel]
        peak_values = []
        for location in peak_location:
            peak_values.append(peak_data[location])
        output_dict[base] = peak_values

    traces = []
    for i in range(len(primary_seq)):
        trace = []
        for base in base_order:
            trace.append(output_dict[base][i])
        # 标准化
        trace = list(map(lambda x: x/sum(trace), trace))
        traces.append(trace)
    return {'Primary_seq':primary_seq,'Quality_scores':quality_scores,'Base_order':base_order,'Peak_values':output_dict,'traces_data':traces}


def find_max_high_quality_window(quality_scores,window_size=30,quality_line=40):
    # 平均质量分数
    window_mean_quality = []
    high_quality_window = []
    is_high_quality = False
    for i in range(len(quality_scores)-(window_size-1)):
        # 每30bp求平均质量分数
        mean_quality = sum(quality_scores[i:i+window_size])/window_size
        window_mean_quality.append(mean_quality)
        if is_high_quality:
            if mean_quality>quality_line:
                pass
            else:
                high_quality_end = i
                is_high_quality = False
                high_quality_window.append([high_quality_start,high_quality_end])

        else:
            if mean_quality>quality_line:
                is_high_quality = True
                high_quality_start = i
    # print(window_mean_quality)
    # print(high_quality_window)

    if len(high_quality_window)==1:
        max_window = high_quality_window[0]
        return max_window

    elif len(high_quality_window)>1:
        max_window = sorted(high_quality_window,key=lambda window:window[1]-window[0],reverse = True)[0]
        return max_window
    else:
        return False


def find_target_in_ctrl(guide_label, guide_seq, ctrl_seq):
    revcomp_guide = DNA_rev_complement(guide_seq)

    if guide_seq in ctrl_seq:
        guide_start = ctrl_seq.index(guide_seq)
        guide_end = guide_start + len(guide_seq)
        # cut_offset = len(guide_seq) - 3
        direction = "fw"
        cutsite = guide_end - 3
        seq = guide_seq
        # check PAMs
        if ctrl_seq[guide_end + 1:guide_end + 3] != 'GG':
            warnings.warn('No PAMs')

    elif revcomp_guide in ctrl_seq:
        guide_start = ctrl_seq.index(revcomp_guide)
        guide_end = guide_start + len(revcomp_guide)
        # cut_offset = 3
        cutsite = guide_start + 3
        direction = "rev"
        seq = revcomp_guide
        # check PAMs
        if ctrl_seq[guide_start - 3:guide_start - 1] != 'CC':
            warnings.warn('No PAMs')
        # print(cutsite)
    else:
        return None

    return {'direction':direction,'cutsite':cutsite,'guide_start':guide_start,'guide_end':guide_end,'guide_seq':seq,'guide_label':guide_label}

def find_alignment_window(high_quality_window,cutsite,indel_max_size):

    aln_window_start = high_quality_window[0]
    aln_window_end = cutsite - indel_max_size - 10

    if aln_window_end not in range(high_quality_window[0],high_quality_window[1]):
        return None

    if aln_window_end - aln_window_start < 40:
        aln_window_start = max(1,aln_window_end - 40)

    return (aln_window_start,aln_window_end)

def wtire_alignment(aln,filename):

    align = Align.MultipleSeqAlignment([])
    align.add_sequence('control', aln[0])
    align.add_sequence('experiment', aln[1])
    align_clustal = align.format("clustal").split('\n', 2)[2]

    with open(filename, 'w') as align_output:
        align_output.write(align_clustal)


def single_guide_genotype(guide,wt_data,type,del_before=0,del_after=0,insert=0):

    wt_primary_seq = list(wt_data['Primary_seq'])
    wt_trace_data = copy(wt_data['traces_data'])
    cutsite = guide['cutsite']
    label = guide['guide_label']

    if type == 'delete':
        delete_bases = list(range(cutsite-del_before,cutsite))+list(range(cutsite,cutsite+del_after))

        dac = 0
        for i in delete_bases:
            wt_primary_seq[i] = '-'
            del wt_trace_data[i-dac]
            dac += 1

        return {'primary_seq':wt_primary_seq, 'sequence':(''.join(wt_primary_seq)).replace('-',''), 'trace_data':wt_trace_data,'cutsite':cutsite, 'change':-len(delete_bases),'type':'{}[{}]'.format(-(len(delete_bases)),label)}

    elif type == 'insert':

        wt_primary_seq = wt_primary_seq[:cutsite] + ['n']*insert + wt_primary_seq[cutsite:]
        wt_trace_data = wt_trace_data[:cutsite] + [[0.25,0.25,0.25,0.25]]*insert + wt_trace_data[cutsite:]

        return {'primary_seq': wt_primary_seq, 'sequence':(''.join(wt_primary_seq)).replace('-',''), 'trace_data': wt_trace_data, 'cutsite': cutsite,'change': insert, 'type': '{}[{}]'.format(+(insert), label)}
        # print(wt_primary_seq,wt_trace_data)

def mutilpe_guide_genotype(guide1,guide2,wt_data,type,cut1_before=0,cut1_after=0,cut2_before=0,cut2_after=0,cut1_insert=0,cut2_insert=0,dropout=False):

    wt_primary_seq = list(wt_data['Primary_seq'])
    wt_trace_data = copy(wt_data['traces_data'])
    cutsite1 = guide1['cutsite']
    label1 = guide1['guide_label']
    cutsite2 = guide2['cutsite']
    label2 = guide2['guide_label']

    if type == 'delete':
        if dropout:
            cut1_bases = list(range(cutsite1 - cut1_before, cutsite1))
            cut2_bases = list(range(cutsite2, cutsite2 + cut2_after))
            delete_bases = cut1_bases + list(range(cutsite1,cutsite2)) + cut2_bases

        else:
            # 不要影响到另一个切割位点
            cut1_bases = list(range(cutsite1 - cut1_before, cutsite1)) + list(range(cutsite1, min(cutsite1 + cut1_after, cutsite2)))
            cut2_bases = list(range(max(cutsite2 - cut2_before, cutsite1 + 1), cutsite2)) + list(range(cutsite2, cutsite2 + cut2_after))
            # 去除重复项（两个切割位点相距太近时会有重复）
            delete_bases = list(set(cut1_bases + cut2_bases))

        dac = 0
        for i in delete_bases:
            wt_primary_seq[i] = '-'
            del wt_trace_data[i - dac]
            dac += 1

        # 删除了多少个碱基，就在序列后面补多少个N，保持序列长度不变
        wt_primary_seq += 'n' * len(delete_bases)
        wt_trace_data += [[0.25,0.25,0.25,0.25]] * len(delete_bases)

        return {'primary_seq': wt_primary_seq, 'sequence':(''.join(wt_primary_seq)).replace('-',''), 'trace_data': wt_trace_data, 'cutsite': (cutsite1,cutsite2),
                'change': -len(delete_bases), 'type': '-{} -{}[{}],-{}[{}]'.format(len(delete_bases),len(cut1_bases), label1, len(cut2_bases), label2)}

    if type == 'insert':
        if dropout:
            # 将切割位点之间的碱基全部敲除
            delete_bases = list(range(cutsite1,cutsite2))
            dac = 0
            for i in delete_bases:
                wt_primary_seq[i] = '-'
                del wt_trace_data[i - dac]
                dac += 1

            # 然后在第一个切割位点切面加N
            wt_primary_seq = wt_primary_seq[:cutsite1] + ['n'] * cut1_insert + wt_primary_seq[cutsite1:]
            wt_trace_data = wt_trace_data[:cutsite1] + [[0.25, 0.25, 0.25, 0.25]] * cut1_insert + wt_trace_data[cutsite1:]

            # 删除了多少个碱基，就在序列后面补多少个N，保持序列长度不变
            wt_primary_seq += 'n' * len(delete_bases)
            wt_trace_data += [[0.25, 0.25, 0.25, 0.25]] * len(delete_bases)

            base_change = cut1_insert+cut2_insert-len(delete_bases)
            summary = '{} +{}[{}] +{}[{}]'.format(base_change, cut1_insert, label1, cut2_insert, label2)


        else:
            wt_primary_seq = wt_primary_seq[:cutsite1] + ['n'] * cut1_insert + wt_primary_seq[cutsite1:cutsite2] + ['n'] * cut2_insert + wt_primary_seq[cutsite2:]
            wt_trace_data = wt_trace_data[:cutsite1] + [[0.25, 0.25, 0.25, 0.25]] * cut1_insert + wt_trace_data[cutsite1:cutsite2] + [[0.25, 0.25, 0.25, 0.25]] * cut2_insert + wt_trace_data[cutsite2:]
            base_change = cut1_insert+cut2_insert
            summary = '{} +{}[{}] +{}[{}]'.format(base_change, cut1_insert, label1, cut2_insert, label2)

        return {'primary_seq': wt_primary_seq, 'sequence': (''.join(wt_primary_seq)).replace('-', ''), 'trace_data': wt_trace_data,
                'cutsite': (cutsite1 + cut1_insert, cutsite2 + cut1_insert + cut2_insert), 'change': base_change, 'type': summary}

def write_all_genotypes(genotypes,indel_max_size,filename):

    with open(filename, 'w') as f:
        for genotype in genotypes:
            if isinstance(genotype['cutsite'], int):
                cutsite = genotype['cutsite']
                genotype['primary_seq'].insert(cutsite, '|')
                start = cutsite - indel_max_size - 5
                end = cutsite + indel_max_size +20
                f.write(''.join(genotype['primary_seq'][start:end])+'  '+genotype['type']+'\n')

            elif isinstance(genotype['cutsite'], tuple):
                cutsite1 = genotype['cutsite'][0]
                cutsite2 = genotype['cutsite'][1]
                genotype['primary_seq'].insert(cutsite2, '|')
                genotype['primary_seq'].insert(cutsite1, '|')
                start = cutsite1 - indel_max_size - 5
                end = cutsite2 + indel_max_size +20
                f.write(''.join(genotype['primary_seq'][start:end]) + '  ' + genotype['type'] + '\n')


