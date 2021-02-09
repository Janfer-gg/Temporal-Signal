import os
import traceback

from Ubigene.ICE.sanger_analysis import SangerAnalysis,ICEResult

# control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK19-004-Npnt/WT.ab1')
# sample_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK19-004-Npnt/PK200908-02-A05-PK593.ab1')
# guide = 'AGGTGCCCTATCGTGTTCCA,TGCAGCTAGCTCCCGATGGG'

# control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK20-055-METTL3/WT.ab1')
# sample_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/CK20-055-METTL3/PK200908-02-A07-PK374.ab1')
# guide = 'GGGCTGTCACTACGGAAGGT,AGCATCAGTGGGCAATGTTA'

# control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/PK1031/PK201114-03-A02-PK1031_WT.ab1')
# sample_path = os.path.abspath('C:/Users/41518/Desktop/ICE/run_exprisement(1)/PK1031/PK210116-06-H07-PK1031.ab1')
# guide = 'CAGCTCCTGTGCGGCCCCCT,aaccagatgaaagagcgcca'

control_path = os.path.abspath('C:/Users/41518/Desktop/ICE/ice/H9c2-Myh7-PK368_WT.ab1')
sample_path = os.path.abspath('C:/Users/41518/Desktop/ICE/ice/PK210107-08-A01-PK368.ab1')
guide = 'GTTATCATTCCGAACTGTCT'
base_outputname = './good_example/'
donor = 'GGCACCTTGGAAGATCAAATCATCCAAGCCAACCCCGCTCTGGAGGCCTTTGGCAATGCCCAAACAGTTCGGAATGATAACTCCTCCCGATTTGTGAGTGATACCCCACCTTGAACTCGGGAC'
# donor = None
verbose = True


def single_sanger_analysis(control_path, sample_path, base_outputname, guide, donor=None, verbose=False,allprops=True):

    # 上传文件异常处理
    if control_path is None or not os.path.exists(control_path):
        raise Exception('Control @ {} not found'.format(control_path))

    if not os.path.exists(sample_path):
        raise Exception('Experiment sample @ {} not found'.format(sample_path))

    # 创建生成结果目录
    base_dir = os.path.join(*os.path.split(os.path.abspath(base_outputname))[:-1])
    if verbose:
        print('Base dir: %s' % base_dir)

    if not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    # 初始化SangerAnalysis类型
    sa = SangerAnalysis(verbose=verbose)
    sa.allprops = allprops

    # 获取属性 control_sample，edited_sample，gRNA_sequences，base_outputname，indel_max_size，donor，allprops
    # control_sample和edited_sample是SangerObject类型，拥有属性data：测序文件的数据，phred_scores：测序质量分数，path：文件路径，basename：文件名
    sa.initialize_with(control_path=control_path,
                       edited_path=sample_path,
                       gRNA_sequences=guide,
                       indel_max_size=15,
                       base_outputname=base_outputname,
                       donor=donor,
                       allprops=allprops)


    try:
        sa.analyze_sample()
        # print(sa.alignment_window)
        return sa.results.to_json(sa.guide_targets, sa.warnings)

    except Exception as e:
        results = ICEResult()
        print('Exception Caught!')
        traceback.print_exc()
        if isinstance(e, KeyError):
            return results.to_json(sa.guide_targets, [';'.join(sa.warnings)])
        else:
            return results.to_json(sa.guide_targets, [str(e)])

if __name__ == '__main__':
    results = single_sanger_analysis(control_path=control_path,
                                 sample_path=sample_path,
                                 base_outputname=base_outputname,
                                 guide=guide,
                                 donor=donor,
                                 verbose=verbose)
    print(results)