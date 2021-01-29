from scipy.stats.stats import pearsonr
from sklearn import linear_model

from Ubigene.ice_analysis.function import get_data,find_target_in_ctrl,wtire_alignment,single_guide_genotype,mutilpe_guide_genotype,find_max_high_quality_window,write_all_genotypes
from Ubigene.ice_analysis.sequence import gRNA_normalize
from Bio import pairwise2
import re
import itertools
import time
import numpy as np
class Analysis:

    def __init__(self,control_path, experiment_path, gRNA_sequences,indel_max_size,outputpath):

        self.control_filepath = control_path  # the file path
        self.control_data = get_data(control_path)

        self.experiment_filepath = experiment_path  # the file path
        self.experiment_data = get_data(experiment_path)

        # the guide sequence
        self.gRNA_sequences = gRNA_normalize(gRNA_sequences)
        self.guide_targets = self.find_targets()   # gRNA切割位点

        self.outputpath = outputpath
        self.indel_max_size = indel_max_size

        self.alignment_window = None  # 比对窗口
        self.all_aligned_seq = None   # 全序列比对后的两条序列
        self.window_aligned_seq = None  # 窗口序列比对后的两条序列
        self.alignment_pairs_index =None   # 碱基对索引

        self.genotypes = None    # 基因型列表
        self.inference_window = None  # 分析窗口
        self.analysis_trace_martix = None
        self.analysis_trace_vector = None

        self.R_squared = None


        self.valid_configuration = False
        self.donor_odn = None
        self.recombination_changed_bases = []

        # alignment of donor with ctrl
        self.donor_alignment = None

        self.alignment_sequence = None

        # if false, the inference matrix will be filled with 0s and 1s based on the sequence only
        self.inference_sequence = None

        self.indel_list = None

        # self.results = ICEResult()

        self.verbose = False
        self.use_partial_for_insertions = True
        self.debug = False
        self.warnings = []
        self.MIN_ALN_WINDOW_SIZE = 50  # minimum length of high quality sequence in sample required for alignment

        self.plot_guide_line = True  # False if we want to use css for denoting where the guide lines are

    def find_targets(self):
        guide_targets = []
        for guide_idx, guide in enumerate(self.gRNA_sequences):
            label = "g{}".format(guide_idx + 1)
            g_t = find_target_in_ctrl(label, guide, self.control_data['Primary_seq'])
            if g_t:
                guide_targets.append(g_t)
            else:
                return None
        # 根据切割位点排序
        guide_targets = sorted(guide_targets, key=lambda x: x['cutsite'])
        return guide_targets

    def all_align(self):
        alignments = pairwise2.align.localms(self.control_data['Primary_seq'],self.experiment_data['Primary_seq'],2, -1, -3, -1)
        aln = alignments[0]
        self.all_aligned_seq = aln
        wtire_alignment(aln,self.outputpath + "all_aligned.txt")


    def window_align(self):

        try:
            alignments = pairwise2.align.localms(self.control_data['Primary_seq'][self.alignment_window[0]:self.alignment_window[1]],
                                                 self.experiment_data['Primary_seq'], 2, -1, -2, -1)
            aln = alignments[0]
        except:
            return 'No alignment'

        self.window_aligned_seq = aln
        wtire_alignment(aln,self.outputpath + "window_aligned.txt")

        # 比对分数
        aln_score = aln[2]
        window_size = self.alignment_window[1]-self.alignment_window[0]
        aln_score_normalized = aln_score / (window_size * 2) * 100

        if aln_score_normalized < 50:
            return "Poor alignment"

        else:
            return 'succeed'


    def make_base_index(self):

        seq1 = self.window_aligned_seq[0]
        seq2 = self.window_aligned_seq[1]

        # print(seq2)

        control_index_first = self.alignment_window[0]
        control_index_last = self.alignment_window[1]

        first_align_index = re.search('[A-Z]',seq1).start()
        last_align_index = re.search('[A-Z\-]+[A-Z]',seq1).end()

        experiment_index =[]
        e_index = 0
        try:
            for i in range(len(seq2)):
                if seq2[i] != '-':
                    experiment_index.append(e_index)
                    e_index += 1
                else:
                    experiment_index.append(None)

            control_index_inside = []
            c_index = control_index_first - 1
            for i in range(first_align_index,last_align_index):
                if seq1[i] != '-':
                    c_index += 1
                    control_index_inside.append(c_index)
                else:
                    control_index_inside.append(None)

            if control_index_first < first_align_index:
                control_index_before = [None]*(first_align_index-control_index_first)
                control_index_before += [i for i in range(control_index_first)]
                control_index_after = [i for i in range(control_index_last, len(self.control_data['Primary_seq']))]
            else:
                control_index_before = [i for i in range(control_index_first-first_align_index,control_index_first)]
                control_index_after = [i for i in range(control_index_last, len(self.control_data['Primary_seq']))]

            control_index = control_index_before + control_index_inside + control_index_after

            alignment_pairs_index = []
            for id1,id2 in zip(control_index,experiment_index):
                alignment_pairs_index.append([id1,id2])

            self.alignment_pairs_index = alignment_pairs_index
            return 'succeed'
        except:
            return False

    def generate_genotype(self):

        genotypes = []
        start = time.time()

        # single cut cases
        for guide in self.guide_targets:
            if guide['cutsite'] - self.indel_max_size < 0 or guide['cutsite'] + self.indel_max_size > len(self.control_data['Primary_seq']):
                print('warning: 敲除范围太大')

            # 敲除
            deletion_type_iterator = itertools.product(range(self.indel_max_size + 1),repeat=2)
            for item in deletion_type_iterator:
                del_before = item[0]
                del_after = item[1]

                genotype_single_delete = single_guide_genotype(guide,self.control_data,type='delete',del_before=del_before,del_after=del_after)
                genotypes.append(genotype_single_delete)
                # print(''.join(genotype_single_delete['primary_seq']))

            # 敲入
            for insertion in range(self.indel_max_size + 1):
                genotype_single_insert = single_guide_genotype(guide,self.control_data,type='insert',insert=insertion)
                genotypes.append(genotype_single_insert)

        # mutiple cut cases
        for dual_guide in itertools.combinations(self.guide_targets,2):
            guide1 = dual_guide[0]
            guide2 = dual_guide[1]
            if guide1['cutsite']>=guide2['cutsite']:
                print("Warning: cutsite1 >= cutsite2")

            deletion_case_iterator = itertools.product(range(5), repeat=4)

            # deletion case
            for item in deletion_case_iterator:

                genotype_dual_delete = mutilpe_guide_genotype(guide1,guide2,self.control_data,type='delete',cut1_before=item[0],
                                                              cut1_after = item[1],cut2_before = item[2],cut2_after = item[3])
                genotypes.append(genotype_dual_delete)

            # dropout deletion case
            deletion_dropout_iterator = itertools.product(range(5), repeat=2)
            for item in deletion_dropout_iterator:

                genotype_dual_delete_dropout = mutilpe_guide_genotype(guide1,guide2,self.control_data,type='delete',cut1_before=item[0],
                                                              cut2_after = item[1],dropout=True)
                genotypes.append(genotype_dual_delete_dropout)

            # insertion case
            insertion_case_iterator = itertools.product(range(3), repeat=2)
            for item in insertion_case_iterator:

                genotype_dual_insert = mutilpe_guide_genotype(guide1,guide2,self.control_data,type='insert',cut1_insert=item[0],
                                                              cut2_insert = item[1],dropout=False)
                genotypes.append(genotype_dual_insert)

            # dropout insertion case
            for item in range(3):
                genotype_dual_insert_dropout = mutilpe_guide_genotype(guide1, guide2, self.control_data, type='insert', cut1_insert=item, dropout=True)

                genotypes.append(genotype_dual_insert_dropout)
        # print(len(genotypes))
        write_all_genotypes(genotypes, self.indel_max_size, self.outputpath + 'all genotypes.txt')
        seen = []
        self.genotypes = list(filter(lambda x: seen.append(x['sequence']) is None if x['sequence'] not in seen else False, genotypes))
        # print(len(self.genotypes))
        # end = time.time()
        # print("Running time: %s seconds" % (end - start))

    def find_inference_window(self):

        inference_window_start = self.alignment_window[1]
        last_guide_cutsite = max(map(lambda x: x['cutsite'], self.guide_targets))
        min_indel_sequence_length = min(map(lambda x: len(x['sequence']), self.genotypes))
        high_quality_window = find_max_high_quality_window(self.control_data['Quality_scores'], window_size=10, quality_line=30)

        if high_quality_window:
            if high_quality_window[1] > last_guide_cutsite:
                inference_window_end = min(high_quality_window[1], last_guide_cutsite+100, min_indel_sequence_length)

            else:
                print("ctrl quality too low, not found high_quailty_inference_window")
                inference_window_end = min(last_guide_cutsite+100, min_indel_sequence_length)
        else:
            print("not found high_quailty_inference_window")
            inference_window_end = min(last_guide_cutsite + 100, min_indel_sequence_length)


        if inference_window_end - last_guide_cutsite < 0:
            print("敲除范围太大，敲除后的序列太短 ")

        else:
            print("inference_window after cutsite is {}bp".format(inference_window_end-last_guide_cutsite))

        self.inference_window=(inference_window_start,inference_window_end)

    # 基因型峰值矩阵
    def inference_trace_martix(self):
        trace_list = []
        for genotype in self.genotypes:
            trace_list.append(list(itertools.chain.from_iterable(genotype['trace_data'][self.inference_window[0]:self.inference_window[1]])))
        trace_martix = np.asarray(trace_list)
        # print(trace_martix)
        self.analysis_trace_martix = trace_martix.T
    # 实验组峰值向量
    def exp_trace_vector(self):
        for i in range(len(self.alignment_pairs_index)):
            if self.alignment_pairs_index[i][0] == self.inference_window[0]:
                window_start = self.alignment_pairs_index[i][1]
            elif self.alignment_pairs_index[i][0] == self.inference_window[1]:
                window_end = self.alignment_pairs_index[i][1]

        trace_list = list(itertools.chain.from_iterable(self.experiment_data['traces_data'][window_start:window_end]))

        self.analysis_trace_vector = np.asarray(trace_list)

    def caluate_coefficient(self):

        try:
            lasso_model = linear_model.Lasso(alpha=0.01, positive=True)
            lasso_model.fit(self.analysis_trace_martix, self.analysis_trace_vector)
            coef_matrix = lasso_model.coef_
            predicted_vector = np.dot(self.analysis_trace_martix, coef_matrix)

        except :
            raise ('Fail to fit lasso model')

        # print(list(predicted))
        # print(list(self.analysis_trace_vector))
        (r, p_value) = pearsonr(predicted_vector, self.analysis_trace_vector)
        # fit_r 皮尔逊相关系数
        r_squared = r ** 2
        print("R_SQUARED {}".format(r_squared))

        coef_total = coef_matrix.sum()
        # here we normalize the relative abundance of each possible indel
        # 每种基因型的相对丰度
        coef_rel = coef_matrix/coef_total * r_squared * 100
        # 向下取整
        coef_rel_fl = np.floor(coef_rel)
        # print(coef_rel_fl)
        total_dif = round(r_squared * 100) - np.sum(coef_rel_fl)

        # 按差值大小排序
        order = np.argsort(coef_rel-coef_rel_fl)[::-1]
        # 补差值
        for i in order:
            coef_rel_fl[i] += 1
            total_dif -= 1
            if total_dif <= 0:
                break

        self.R_squared = round(r_squared * 100)
        for genotype,coef in zip(self.genotypes,coef_rel_fl):
            genotype['coef'] = coef

        # 排序
        self.contributions = sorted(self.genotypes,key=lambda x:x['coef'],reverse=True)
        # 去掉coef=0
        self.contributions=list(filter(lambda x:x['coef']!=0,self.contributions))


    def write_contributions(self,filename):

        with open(filename, 'w') as f:
            for contrib in self.contributions:
                if isinstance(contrib['cutsite'], int):
                    cutsite = contrib['cutsite']
                    contrib['primary_seq'].insert(cutsite, '|')
                    start = cutsite - self.indel_max_size - 5
                    end = cutsite + self.indel_max_size + 20
                    f.write(str(contrib['coef']) + ' \t ' + str(contrib['type']) + ' \t ' +  ''.join(contrib['primary_seq'][start: end]) + '\n')

                elif isinstance(contrib['cutsite'], tuple):
                    cutsite1 = contrib['cutsite'][0]
                    cutsite2 = contrib['cutsite'][1]
                    contrib['primary_seq'].insert(cutsite2, '|')
                    contrib['primary_seq'].insert(cutsite1, '|')
                    start = cutsite1 - self.indel_max_size - 5
                    end = cutsite2 + self.indel_max_size + 20
                    f.write(str(contrib['coef']) + ' \t ' + str(contrib['type']) + ' \t ' + ''.join(contrib['primary_seq'][start:end]) + '\n')

