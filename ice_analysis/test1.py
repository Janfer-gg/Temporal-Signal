import os
from Ubigene.ice_analysis.function import input_check,find_max_high_quality_window,find_alignment_window
from Ubigene.ice_analysis.SangerObject import Analysis

# control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK19-004-Npnt/WT.ab1')
# experiment_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK19-004-Npnt/PK200908-02-A05-PK593.ab1')
# gRNA = ['AGGTGCCCTATCGTGTTCCA', 'TGCAGCTAGCTCCCGATGGG']

# control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK20-055-METTL3/WT.ab1')
# experiment_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK20-055-METTL3/PK200908-02-A07-PK374.ab1')
# gRNA = ['GGGCTGTCACTACGGAAGGT', 'AGCATCAGTGGGCAATGTTA']

control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/ice/H9c2-Myh7-PK368_WT.ab1')
experiment_path = os.path.abspath('C:/Users/41518/Desktop/ICE/ice/PK210107-08-A01-PK368.ab1')
gRNA = ['GTTATCATTCCGAACTGTCT']
outputpath = 'good_example/'
donor = 'GGCACCTTGGAAGATCAAATCATCCAAGCCAACCCCGCTCTGGAGGCCTTTGGCAATGCCCAAACAGTTCGGAATGATAACTCCTCCCGATTTGTGAGTGATACCCCACCTTGAACTCGGGAC'
# donor = None
verbose = True
indel_max_size = 15


# 主函数
def single_analysis(control_path, experiment_path, outputpath, gRNA, indel_max_size, donor=None, verbose=False, allprops=True):

    # 输入处理
    check_logic,check_result = input_check(control_path, experiment_path, outputpath, gRNA)
    if check_logic == False:
        return check_result

    sa = Analysis(control_path,experiment_path,gRNA,indel_max_size,outputpath,donor)

    # 质量检测（实验组）
    quality_check_experiment = find_max_high_quality_window(sa.experiment_data['Quality_scores'])
    if not quality_check_experiment:
        return 'experiment sample quality scores too low'

    # 找切割位点
    if not sa.guide_targets:
        return 'gRNA not found in control sample'

    if donor:
        sa.recombination()

    # 质量检测（对照组）
    quality_check_control = find_max_high_quality_window(sa.control_data['Quality_scores'],quality_line=30)
    if not quality_check_control:
        return 'control sample quality scores too low'
    else:
        # 比对窗口
        alignment_window = find_alignment_window(quality_check_control,sa.guide_targets[0]['cutsite'],indel_max_size)

        if not alignment_window:
            return 'control sample quality scores too low, not found alignment window'
        else:
            if sa.donor.donor_seq:

                sa.alignment_window = (alignment_window[0],min(sa.donor.donor_pos[0],alignment_window[1]))
            else:
                sa.alignment_window = alignment_window
    # 全局比对，输出文件
    sa.all_align()
    # 局部比对，输出文件
    align_reuslt = sa.window_align()
    if align_reuslt != 'succeed':
        return align_reuslt
    # 建立索引
    make_index_result = sa.make_base_index()
    if not make_index_result:
        return 'Fail to make index from control to experiment'

    # 构建基因型
    sa.generate_genotype()
    # write_genotypes(sa.genotypes,indel_max_size,outputpath+'all genotypes.txt')
    print("analyzing {} number of edit proposals".format(len(sa.genotypes)))
    # [print(str(ind['change']) + '__' + str(len(ind['sequence']))) for ind in sa.genotypes]

    sa.find_inference_window()
    sa.inference_trace_martix()
    sa.exp_trace_vector()
    sa.calculate_coefficient()

    sa.write_contributions(outputpath + 'contrib.txt')
    sa.write_contributions_json()
    sa.write_result()
    return 'succeed'

if __name__ == '__main__':
    result = single_analysis(control_path, experiment_path, outputpath, gRNA, indel_max_size,donor, verbose)
    print(result)
