import os
from itertools import combinations
from collections import defaultdict
import warnings
from scipy.optimize import nnls
from scipy.stats.stats import pearsonr
from sklearn import linear_model
import re
import numpy as np
from Bio import SeqIO
from Ubigene.ICE.misc import running_mean, find_regions_greater_than

from Ubigene.ICE.ice.classes.edit_proposal_creator import EditProposalCreator
from Ubigene.ICE.ice.classes.guide_target import GuideTarget
from Ubigene.ICE.ice.classes.pair_alignment import PairAlignment
from Ubigene.ICE.ice.classes.proposal_base import ProposalBase
from Ubigene.ICE.ice.outputs.create_discordance_indel_files import generate_discordance_indel_files
from ice.outputs.create_json import write_individual_contribs, write_contribs_json, write_all_proposals_json
from ice.outputs.create_trace_files import generate_trace_files

from Ubigene.ICE.sequence import RNA2DNA, reverse_complement


class SangerObject:
    """

    """

    PRIMARY_BASE_CALL = 'PBAS1'
    BASE_ORDER = 'FWO_1'
    PEAK_LOCATIONS = 'PLOC1'
    CHANNELS = ['DATA9', 'DATA10', 'DATA11', 'DATA12']

    def __init__(self, verbose=False):
        self.data = None
        # self.left_trim_bases = 5
        # self.right_trim_bases = 5
        # self.alignment_offset = 0
        self.external_base_calls = None
        self.path = None
        self.basename = None

    def estimate_quality_scores(self):
        # 一种质量评估方法
        '''
        a good enough way to estimate quality scores
        As ICE uses the quality scores for determining windows for alignment, we just need a rough estimate
        of base quality.  If the primary color is 90% of the signal, the phred score is 50% of the maximum.
        If the primary color is 80% of the total signal at that position, the phred score is 0% of the maximum.
        :return: list of phred scores
        '''

        MAX_PHRED = 60
        peak_values = self.get_peak_values()
        phred_scores = []
        for idx, base in enumerate(peak_values['C']):
            sum_values = 0
            max_value = 0
            for color in 'ACGT':
                sum_values += peak_values[color][idx]
                if peak_values[color][idx] > max_value:
                    max_value = peak_values[color][idx]

            if sum_values < 100 or sum_values > 5000:
                phred_score = 0
            else:
                max_value_percent = max_value / sum_values * 100

                # linear relationship of PHRED score 0 to 60 and max_value_percent 80 to 100
                phred_score = max(3 * max_value_percent - 240, 0)
            # check if absolute intensities are reasonable

            phred_scores.append(phred_score)
        return phred_scores

    def initialize_from_path(self, ab1_file_path, *args, **kwargs):
        '''
        :param ab1_file_path:
        :param args:
        :param kwargs:
        :return:

        'FWO_1': 'GATC' contains the base order
        DATA9: G, DATA10: A, DATA11: T, DATA12: C
        good documentation of what all the data channels are here:
        https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py
        '''
        self.path = ab1_file_path
        self.basename = os.path.basename(ab1_file_path)
        record = SeqIO.read(ab1_file_path, 'abi')

        # In Biopython 1.74, some strings are converted to bytes in py3
        # Convert strings with bytes-object back to regular string
        # 解析测序数据
        traces_ = record.annotations['abif_raw']
        traces = {}
        for k, v in traces_.items():
            if isinstance(v, bytes):
                v = v.decode()
            else:
                v = v
            traces[k] = v
        self.data = traces

        # 质量分数
        phred_scores = []
        for c in traces['PCON2']:
            phred_scores.append(ord(c))
        self.phred_scores = phred_scores
        override_quality_score = kwargs.get('override_quality_score', False)
        # 如果质量分数全都是0，则自己评估质量分数
        if np.count_nonzero(phred_scores) == 0 or override_quality_score:
            self.phred_scores = self.estimate_quality_scores()

    def find_alignable_window(self, window_size=30, QUAL_CUTOFF=40):
        '''
        This code finds the optimal window for aligning the edited trace to the reference sequence
        1) It looks at the phred scores and calculates a moving average
        2) It then selects the largest window where the moving average is greater than QUAL_CUTOFF
        3) This segment and other data is then returned in a dictionary
        all coordinates are absolute with respect to the sanger object
        '''

        # TODO: are there any edge cases where an edit is so good that the alignable window would extend into
        # regions where there is mismatch between sample and control?

        phreds = np.asarray(self.phred_scores)
        # find a running mean
        window_size = window_size
        # 每30bp求平均质量分数
        windowed_phred = running_mean(phreds, window_size)
        # print(windowed_phred)
        # find regions with windowed mean > QUAL_CUTOFF
        # 找到平均质量分数大于40分的区域（一个或多个）
        y = find_regions_greater_than(windowed_phred, QUAL_CUTOFF)
        # print(y)
        max_window_size = 0
        max_window = None
        for region in y:
            w_size = region[1] - region[0]
            if w_size > max_window_size:
                max_window_size = w_size
                max_window = (region[0], region[1])

        return {
            'regions': y,
            'max_window': max_window,
            'windowed_phred': windowed_phred
        }

    @property
    def primary_base_calls(self):
        return self.data[SangerObject.PRIMARY_BASE_CALL]

    @property
    def aligned_bases(self):
        return self.primary_base_calls

    @property
    def base_order(self):
        return self.data[SangerObject.BASE_ORDER]

    @property
    def peak_locations(self):
        """ Returns the peak locations for the final sequence

        :return:
        """
        return list(self.data[SangerObject.PEAK_LOCATIONS])

    # 峰所在位置的GACT各自的峰值
    def get_peak_values(self):
        output_dict = {}
        for n, channel in enumerate(SangerObject.CHANNELS):
            peak_data = np.array(self.data[channel])
            # peak_values = peak_data[self.peak_locations]
            peak_values = peak_data[list(self.data[SangerObject.PEAK_LOCATIONS])]
            output_dict[self.base_order[n]] = peak_values
        # print(output_dict)
        return output_dict

    def get_trough_position(self, loc1, side="left"):
        if side == "left":
            loc2 = loc1 - 1
        else:
            loc2 = loc1 + 1
        return (self.peak_locations[loc1]
                + self.peak_locations[loc2]) / 2

    @classmethod
    def trace_to_base_calls(cls, trace, base_order):
        BASE_CALL_CUTOFF = 0.75
        assert len(trace) % 4 == 0
        seq_len = int(len(trace) / 4)
        seq = ""
        for base_idx in range(seq_len):
            slice_values = trace[base_idx * 4: base_idx * 4 + 4]
            slice_max = max(slice_values)
            if slice_max > BASE_CALL_CUTOFF:
                max_idx = slice_values.index(slice_max)
                base = base_order[max_idx]
                seq += base
            else:
                seq += 'N'
        return seq


class ICEResult:

    def __init__(self):
        self.control_discordances = []
        self.edited_discordances = []
        self.mean_discord_before = None
        self.mean_discord_after = None
        self.aggregate = None
        self.contribs = None
        self.edit_eff = None
        self.ice_d = None
        self.r_squared = None
        self.hdr_percentage = None
        self.ko_score = None

    def to_json(self, guide_targets, warnings):

        results = {}
        if self.edit_eff is not None:
            results['ice'] = int(round(self.edit_eff * 100, 0))
        else:
            results['ice'] = None
        if self.mean_discord_before is not None:
            results['mean_discord_before'] = round(self.mean_discord_before, 3)
        else:
            results['mean_discord_before'] = None
        if self.mean_discord_after is not None:
            results['mean_discord_after'] = round(self.mean_discord_after, 3)
        else:
            results['mean_discord_after'] = None
        if self.ice_d is not None:
            results['ice_d'] = int(round(self.ice_d, 0))
        else:
            results['ice_d'] = None
        if self.r_squared is not None:
            results['rsq'] = round(self.r_squared, 2)
        else:
            results['rsq'] = None

        results['hdr_pct'] = self.hdr_percentage
        results['ko_score'] = self.ko_score
        results['notes'] = "; ".join(warnings)
        results['guides'] = [vars(g) for g in guide_targets]

        results['status'] = 'succeeded'

        return results


def round_percent(orig_array, r_squared):
    # Scale the array by 100 in order to work with np.floor
    scaled_array = np.array([x * 100 for x in orig_array])

    # floored array
    # 向下取整数，去掉小于1的
    fl_dist = np.floor(scaled_array)
    # Calculate the total difference between the floored array and the r2, this difference must be added back
    # 原本fl_dist的总和等于r2，去掉小于1的之后，求出这个差值
    total_lost = round(r_squared * 100) - np.sum(fl_dist)

    # Rank each value by how much it lost when it was floored
    # 按照损失的值由大到小排序（由于向下取整而损失的值）
    order = sorted([(x, ix) for ix, x in enumerate(scaled_array - np.floor(scaled_array))], reverse=True)
    # Add back what was lost in total, accoridng to how much each value lost.
    # 按顺序补差值，每次+1
    for l in range(int(total_lost)):
        fl_dist[order[l][1]] += 1

    return fl_dist / 100  # set it back to a percentage


class SangerAnalysis:
    """
    1. Align sgRNA guide sequence to control sequence (or to the base calls from the control chromatogram)
    2. Align the upstream portion of the control sequence to the upstream portion of the experimental sample
    3. For each base call in experimental sample, extract relative abundance of each nucleotide
    4. Sum relative abundances of incorrect bases, normalize to 75 -> 100%

    """

    iced_correction_factor = 1.41  # ICE-D regression coefficient
    HDR_OVERLAP_FILTER_CUTOFF = 1  # abs val cutoff to not skip proposals even if they overlap in size w/ HDR proposal

    def __init__(self, verbose=False):

        self.control_chromatogram_path = None  # the file path
        self.control_sample = None  # the Sanger Object

        self.edited_chromatogram_path = None  # the file path
        self.edited_sample = None  # the Sanger Object

        # self.gRNA_sequence = None  # the guide sequence
        self.gRNA_sequences = []
        self.guide_targets = []

        self.valid_configuration = False
        self.donor_odn = None
        self.recombination_changed_bases = []

        # list of EditProposals
        self.proposals = None

        # alignment of donor with ctrl
        self.donor_alignment = None

        self.alignment_window = None  # these coordinates refer to the ctrl
        self.inference_window = None  # these coordinates refer to the ctrl, not the sample
        self.alignment_sequence = None
        self.indel_max_size = None

        self.base_outputname = None

        self.use_ctrl_trace = True  # if true, this will use the peak values from the ctrl trace for inference
        # if false, the inference matrix will be filled with 0s and 1s based on the sequence only
        self.inference_sequence = None
        self.indel_list = None
        self.output_vec = None
        self.results = ICEResult()

        self.r_squared_correction = True

        self.verbose = verbose
        self.use_partial_for_insertions = True
        self.debug = False
        self.warnings = []
        self.MIN_ALN_WINDOW_SIZE = 50  # minimum length of high quality sequence in sample required for alignment

        self.plot_guide_line = True  # False if we want to use css for denoting where the guide lines are
        self.allprops = False  # print out all proposals in the json, even if they have 0 contribution

    def initialize_with(self, control_path, edited_path, gRNA_sequences,
                        alignment_window=None,
                        inference_window=None,
                        indel_max_size=5,
                        base_outputname=None,
                        donor=None,
                        allprops=False):
        """
           indel_size_range=(1,10)

        :param control_path:
        :param edited_path:
        :param gRNA_sequences:
        :param args:
        :param kwargs:
        :return:
        """

        if os.path.exists(control_path):
            # 创建初始化对象
            control_sample = SangerObject()
            control_sample.initialize_from_path(control_path)
            self.control_sample = control_sample
        else:
            raise Exception('Control path does not exist')

        if os.path.exists(edited_path):
            # 创建初始化对象
            edited_sample = SangerObject()
            edited_sample.initialize_from_path(edited_path)
            self.edited_sample = edited_sample
        else:
            raise Exception('Edited path does not exist')

        # gRNA序列不为None且类型为str
        if gRNA_sequences is not None and (isinstance(gRNA_sequences, str)):
            # 将多条序列按逗号分开
            gRNA_sequences = gRNA_sequences.split(",")
            gRNA_sequences = list(map(str.strip, gRNA_sequences))
            gRNA_sequences = list(map(RNA2DNA, gRNA_sequences))  # RNA转DNA
            gRNA_sequences = list(map(lambda x: x.upper(), gRNA_sequences))  # 小写字母转大写字母
            self.gRNA_sequences = gRNA_sequences

        # 没有输入gRNA
        else:
            raise Exception('{} contains invalid gRNA sequences'.format(gRNA_sequences))

        if base_outputname is not '':
            self.base_outputname = base_outputname + "."

        self.inference_window = inference_window
        self.indel_max_size = indel_max_size
        self.donor_odn = donor
        self.allprops = allprops
        return True

    @property
    def multiplex(self):
        """
        :return: bool indicating if sanger analysis is multiplex
        """
        return len(self.guide_targets) > 1

    def quality_check(self):
        # SangerObject类的方法find_alignable_window
        aw = self.edited_sample.find_alignable_window()
        # aw.regions：所有大于40分的区域；aw.max_window：大于40分的区域里最大的一个；aw.windowed_phred：平均质量分数
        if aw['max_window'] is None:
            raise Exception("Sample ab1 %s quality scores too low" % self.edited_sample.basename)

    def sample_cutsite(self, guide_target, alt=False):
        ctrl_cutsite = guide_target.cutsite
        sample_cutsite = self.alignment.ctrl2sample_coords(ctrl_cutsite)

        return sample_cutsite

    # 在控制组序列上找gRNA靶位点
    def find_guide_in_ctrl(self, guide_label, guide_seq, ctrl_seq):
        revcomp_guide = reverse_complement(guide_seq)
        found_seq = None
        if guide_seq in ctrl_seq:
            cut_offset = len(guide_seq) - 3
            orientation = "fwd"
            found_seq = guide_seq
        elif revcomp_guide in ctrl_seq:
            cut_offset = 3
            orientation = "rev"
            found_seq = revcomp_guide
        else:
            raise Exception("guide %s not found in control sequence" % guide_seq)

        if found_seq:
            guide_start = ctrl_seq.index(found_seq)
            guide_end = guide_start + len(found_seq)
            cutsite = guide_start + cut_offset
            ### checking for PAMs
            if orientation == "rev" and ctrl_seq[guide_start - 3:guide_start - 1] != 'CC':
                self.warnings.append("No PAM upstream of guide {}".format(guide_label))

            elif orientation == "fwd" and ctrl_seq[guide_end + 1:guide_end + 3] != 'GG':
                self.warnings.append("No PAM downstream of guide {}".format(guide_label))

            return GuideTarget(
                orientation=orientation,
                cut_offset=cut_offset,
                cutsite=cutsite,
                guide_start=guide_start,
                guide_end=guide_end,
                sequence=guide_seq,
                label=guide_label
            )

    def find_targets(self):
        '''
        This function verifies that the sgRNA sequence is in the control sample or reference,
        it also checks if the reverse complement is in the sequence
        :return:
        '''
        if self.control_sample and self.edited_sample and self.gRNA_sequences:
            for guide_idx, guide in enumerate(self.gRNA_sequences):
                label = "g{}".format(guide_idx + 1)
                g_t = self.find_guide_in_ctrl(label, guide, self.control_sample.primary_base_calls)
                self.guide_targets.append(g_t)
            # 根据切割位点排序
            self.guide_targets = sorted(self.guide_targets, key=lambda x: x.cutsite)

    def check_recombination(self):
        epc = EditProposalCreator(self.control_sample.primary_base_calls,
                                  use_ctrl_trace=True,
                                  sanger_object=self.control_sample)

        # donor太长，警告
        if len(self.donor_odn) > len(self.control_sample.primary_base_calls) * 0.75:
            self.warnings.append(
                "Large Donor of {} bp compared to control sequence of {} bp".format(len(self.donor_odn), len(
                    self.control_sample.primary_base_calls)))

            # int((len(re.split("(-+)", alignment[0])) - 3) / 2)

        try:
            cutsite = self.guide_targets[0].cutsite
            hr_proposal, changed_bases, odn_start_pos, aln = epc.homologous_recombination_proposal(
                cutsite=cutsite, donor_sequence=self.donor_odn)

            self.recombined_seq = hr_proposal.sequence
            self.recombination_changed_bases = changed_bases
            self.odn_start_pos = odn_start_pos
            self.donor_alignment = aln

            n_splits = int((len(re.split("(-+)", aln.all_aligned_seqs[0])) - 3) / 2)

            # 有多个重组位点
            if n_splits > 0:
                self.warnings.append("{} donor integration sites".format(n_splits + 1))


        except Exception as e:
            self.warnings.append("Could not analyze homologous recombination case: {}".format(e))
            self.donor_odn = None

    def find_alignment_window(self):
        """
        Determines the alignment window

        """
        # 平均质量分数大于30的区域
        ctrl_aw = self.control_sample.find_alignable_window(QUAL_CUTOFF=30)
        # print(ctrl_aw)
        if ctrl_aw['max_window'] is None:
            # no quality windows found
            raise Exception("Control ab1 trace quality scores too low")
        start_of_ctrl_aln_window = ctrl_aw['max_window'][0]
        # alignment window should stop upstream of any potential indels
        end_of_ctrl_aln_window = self.guide_targets[0].cutsite - self.indel_max_size - 10

        if end_of_ctrl_aln_window > ctrl_aw['max_window'][1]:
            raise Exception("Control ab1 trace quality too upstream of guide too low")
        if ctrl_aw['max_window'][0] > end_of_ctrl_aln_window:
            raise Exception("Control ab1 %s quality before cutsite is too low" % self.control_sample.basename)
        if end_of_ctrl_aln_window - ctrl_aw['max_window'][0] < 40:
            self.warnings.append("Padding alignment window to be at least 40bp in length or to start of seq")
            start_of_ctrl_aln_window = max(1, end_of_ctrl_aln_window - 40)

        # 暂时忽略
        if self.donor_odn is not None:
            if self.odn_start_pos < end_of_ctrl_aln_window:
                end_of_ctrl_aln_window = self.odn_start_pos

        self.alignment_window = (start_of_ctrl_aln_window, end_of_ctrl_aln_window)

    def calculate_discordance(self):

        # get the peak values at every base
        control_peak_values = self.control_sample.get_peak_values()
        edited_peak_values = self.edited_sample.get_peak_values()

        # ... get the primary sequence; note that the aligned_bases here refers to correctly called bases in
        # the chromatogram
        ctrl_sequence = self.control_sample.aligned_bases
        shortest_sequence = min(len(ctrl_sequence), len(edited_peak_values['T']))

        control_discord_abundances = []
        edited_discord_abundances = []

        for aligned_bases in self.alignment.alignment_pairs:
            ctrl_idx = aligned_bases[0]
            sample_idx = aligned_bases[1]
            if ctrl_idx is None or sample_idx is None or ctrl_idx > shortest_sequence:
                # gap in alignment; gap is in ctrl sample
                continue

            bases = ['A', 'C', 'G', 'T']
            ctrl_base = self.control_sample.primary_base_calls[ctrl_idx]
            # sample_base = self.edited_sample.primary_base_calls[sample_idx]

            if ctrl_base not in bases:

                control_discord_abundances.append(1.0)
                edited_discord_abundances.append(1.0)

            else:
                # calculate the total intensity for this base call.
                total_edited_signal = 0
                total_control_signal = 0

                for ab in bases:
                    total_edited_signal += edited_peak_values[ab][sample_idx]
                    total_control_signal += control_peak_values[ab][ctrl_idx]

                # get the correct signal
                edited_correct_signal = edited_peak_values[ctrl_base][sample_idx]
                control_correct_signal = control_peak_values[ctrl_base][ctrl_idx]

                # now figure out the incorrect signal
                edited_discord_signal = (total_edited_signal - edited_correct_signal)
                control_discord_signal = (total_control_signal - control_correct_signal)

                # and get the relative amount
                edited_discord_signal_relative = (edited_discord_signal / (1.0 * total_edited_signal))
                control_correct_signal = (control_discord_signal / (1.0 * total_control_signal))

                # ta da
                edited_discord_abundances.append(edited_discord_signal_relative)
                control_discord_abundances.append(control_correct_signal)

        final_sequence = ctrl_sequence[0:shortest_sequence]

        guide_location = self.guide_targets[0].guide_start

        cutsite = self.sample_cutsite(self.guide_targets[0], alt="downstream")

        start = self.alignment_window[0]
        end = self.alignment_window[1]
        sample_start = self.alignment.ctrl2sample_coords(start)
        sample_end = self.alignment.ctrl2sample_coords(end)
        before_dsb = edited_discord_abundances[sample_start:sample_end]
        after_dsb = edited_discord_abundances[cutsite:]

        if self.debug:
            print("Cutsite for discord calcs: {}".format(cutsite))

        if not before_dsb or not after_dsb:
            # TODO find more precise description of this error
            self.warnings.append("Unable to compute discordance plots")

        mean_before = np.mean(before_dsb)
        mean_after = np.mean(after_dsb)

        validation_msg = 'discord (aln window): %0.2f after cutsite: %0.2f' % (mean_before, mean_after)

        print(validation_msg)

        self.results.control_discordances = control_discord_abundances
        self.results.edited_discordances = edited_discord_abundances
        self.results.mean_discord_after = mean_after
        self.results.mean_discord_before = mean_before

    def _calculate_inference_window(self):

        # 比对窗口的最右端
        left_offset = self.alignment_window[1]
        # 最靠后的gRNA的切割位点
        last_guide = max(self.guide_targets, key=lambda x: x.cutsite)
        cutsite = last_guide.cutsite
        min_indel_sequence_length = 10000

        # 最靠后的gRNA的切割位点再往后100bp
        MAX_BASES_AFTER_CUTSITE = 100

        # 在proposals的基因型里找到最小长度
        # find minimum length of all generated sequences
        for ind in self.proposals:
            if len(ind.sequence) < min_indel_sequence_length:
                min_indel_sequence_length = len(ind.sequence)
                # min_indel = ind.sequence
        # [print(str(ind.bases_changed) + '__' + str(len(ind.sequence))) for ind in self.proposals]
        # import pdb;pdb.set_trace()
        ctrl_quality_windows = self.control_sample.find_alignable_window(window_size=10, QUAL_CUTOFF=35)

        # we should not be doing any calculations with data from low quality regions
        # we can cutoff the inference window on the right based on quality
        iw_right_boundary = None
        if 'max_window' in ctrl_quality_windows:
            ctrl_aw = ctrl_quality_windows['max_window']

            if ctrl_aw[1] > cutsite:
                iw_right_boundary = min(cutsite + MAX_BASES_AFTER_CUTSITE, ctrl_aw[1],
                                        min_indel_sequence_length)
            else:
                iw_right_boundary = min(cutsite + MAX_BASES_AFTER_CUTSITE,
                                        min_indel_sequence_length)
                self.warnings.append(
                    "Low quality control trace; using cutsite + {} for inference. Results may be poor.".format(
                        MAX_BASES_AFTER_CUTSITE))
        if iw_right_boundary is None:
            self.warnings.append(
                "Low quality control trace; using cutsite + {} for inference. Results may be poor.".format(
                    MAX_BASES_AFTER_CUTSITE))
            iw_right_boundary = cutsite + MAX_BASES_AFTER_CUTSITE

        # as an extra precaution, we want to exclude the last aligned 10bp
        # 为什么是这个？
        last_aligned_ctrl_base = self.alignment.alignment_pairs[self.alignment.last_aligned_pair_idx][0]
        iw_right_boundary = min(iw_right_boundary, last_aligned_ctrl_base - 10)

        inf_len_after_cutsite = iw_right_boundary - cutsite
        # 窗口右端和切割位点的距离小于3倍indel_max_size，警告
        if inf_len_after_cutsite < self.indel_max_size * 3:
            self.warnings.append(
                "Inf. window after cutsite is only {} in length and is less than 3x indel_max_size of {}".format(
                    inf_len_after_cutsite, self.indel_max_size))
        # set inference_window
        self.inference_window = (left_offset, iw_right_boundary)

    def _should_skip_proposal(self, indel_size):
        """
        Indels larger than min cutoff that overlap in size with donor are not considered as an edit proposal
        """
        # Only consider skipping proposals if there is a donor
        if not self.donor_odn:
            return False

        # Don't filter out small indels because they are still likely to occur with HDR
        if abs(indel_size) <= self.HDR_OVERLAP_FILTER_CUTOFF:
            return False

        # Skip proposal if same size as donor
        return indel_size == self.donor_alignment.hdr_indel_size

    def _generate_edit_proposals(self):
        import time
        start = time.time()
        epc = EditProposalCreator(self.control_sample.primary_base_calls,
                                  use_ctrl_trace=True,
                                  sanger_object=self.control_sample)

        proposals = []
        # 暂时忽略
        if self.donor_odn is not None:
            cutsite = self.guide_targets[0].cutsite

            hr_proposal, changed_bases, odn_start_pos, aln = epc.homologous_recombination_proposal(
                cutsite=cutsite, donor_sequence=self.donor_odn)
            proposals.append(hr_proposal)

        # single cut cases
        # 每条gRNA都运行一次

        for guide_target in self.guide_targets:
            cutsite = guide_target.cutsite
            # deletion case
            # [0-indel_max_size]
            deletion_befores = list(range(self.indel_max_size + 1))
            deletion_afters = list(range(self.indel_max_size + 1))

            for deletion_before in deletion_befores:
                for deletion_after in deletion_afters:
                    indel_size = -(deletion_before + deletion_after)

                    # 暂时忽略(donor)
                    if self._should_skip_proposal(indel_size):
                        continue

                    # 返回EditProposal类型
                    # 441次：len(deletion_befores)*len(deletion_afters),其中有一个野生型（deletion_before=deletion_after=0）
                    ep = epc.single_cut_edit_proposal(cutsite, guide_target.label,
                                                      del_before=deletion_before, del_after=deletion_after)
                    proposals.append(ep)
                    # print(ep.bases_changed)
            # [0-indel_max_size]
            insertions = list(range(self.indel_max_size + 1))
            for insertion in insertions:
                # 暂时忽略(donor)
                if self._should_skip_proposal(insertion):
                    continue
                # 20次：len(insertions)
                ep = epc.single_cut_edit_proposal(cutsite, guide_target.label, insertion=insertion)
                proposals.append(ep)
        # multiplex case
        # we limit the deletion sizes here
        deletion_befores = list(range(5))
        deletion_afters = list(range(5))

        insertions = list(range(3))

        # 两两结合(如果只有一条gRNA，不会进入循环)
        for combo in combinations(self.guide_targets, 2):
            guide1 = min(combo, key=lambda x: x.cutsite)
            guide2 = max(combo, key=lambda x: x.cutsite)
            cutsite1 = guide1.cutsite
            cutsite2 = guide2.cutsite
            # 防止两条gRNA相同
            if cutsite1 == cutsite2:
                print("Warning: cutsite1 == cutsite2")
                continue
            label1 = guide1.label
            label2 = guide2.label

            # 5*5*5*5=625
            # 当cut1_del=(0,0)或cut2_del=(0,0)时，直接返回None，其余的都进入第一种情况（576）
            for cut1_before in deletion_befores:
                for cut1_after in deletion_afters:
                    for cut2_before in deletion_befores:
                        for cut2_after in deletion_afters:
                            independent_cut = epc.multiplex_proposal(cutsite1, cutsite2, label1, label2,
                                                                     cut1_del=(cut1_before, cut1_after),
                                                                     cut2_del=(cut2_before, cut2_after))
                            # 625-25-24=576
                            if independent_cut:
                                proposals.append(independent_cut)

            # 移除切割位点之间的全部碱基
            # dropout case
            # cut1_before=0或cut2_after=0时，进入第三种情况（9），当cut1_before和cut2_after都不为0时，进入第一种情况（16）
            for cut1_before in range(5):
                for cut2_after in range(5):
                    dropout = epc.multiplex_proposal(cutsite1, cutsite2, label1, label2, cut1_del=(cut1_before, 0),
                                                     cut2_del=(0, cut2_after), dropout=True)
                    if dropout:
                        proposals.append(dropout)

            # 3*3-1=9种
            # insertion1和insertion2不全为0时，进入第二种情况（8）
            for insertion1 in insertions:
                for insertion2 in insertions:
                    cut_and_insert = epc.multiplex_proposal(cutsite1, cutsite2, label1, label2, cut1_ins=insertion1,
                                                            cut2_ins=insertion2)
                    if cut_and_insert:
                        proposals.append(cut_and_insert)

            # 移除切割位点之间的全部碱基
            # dropout insertion case
            # 当cut1_ins不为0时，进入第二种情况（2）；当cut1_ins为0时，进入第三种情况（1）
            for insertion in insertions:
                dropout_and_insert = epc.multiplex_proposal(cutsite1, cutsite2, label1, label2, cut1_ins=insertion,
                                                            dropout=True)
                if dropout_and_insert:
                    proposals.append(dropout_and_insert)

        end = time.time()
        print("Running time: %s seconds" % (end - start))
        # removing degenerate proposals
        # 过滤去重复（去掉"-"后的序列是否重复）
        # print(len(proposals))
        seen = []
        self.proposals = list(
            filter(lambda x: seen.append(x.sequence) is None if x.sequence not in seen else False, proposals)
        )

    def _generate_coefficient_matrix(self):
        # 所有可能的基因型的数量
        num_proposals = len(self.proposals)

        # inference_window窗口大小
        iw_length = self.inference_window[1] - self.inference_window[0]

        # 零矩阵，行数为num_proposals，列数为4*iw_length
        output_matrix = np.zeros((num_proposals, 4 * iw_length))

        # import pdb; pdb.set_trace()
        # edit_proposal_idx是下标，ep是proposals列表里的元素，每一个ep都是一个基因型
        for edit_proposal_idx, ep in enumerate(self.proposals):
            # base_index是窗口内的下标。（334，553）
            # 遍历碱基
            for base_index in range(self.inference_window[0], self.inference_window[1]):  # （334-552）
                # seq_index是窗口内的相对下标（0-218）
                seq_index = base_index - self.inference_window[0]
                # 获取峰值数据（乘以4是因为trace_data中每4个分别是同一个碱基的GACT峰值）
                # 遍历某个碱基的G、A、C、T
                for color_index in range(4):
                    # 比如想获取第334个碱基的峰值，从334*4的位置开始连续4个数值就是第334个碱基的GACT峰值
                    base_color_index = base_index * 4 + color_index
                    # 同理，output_matrix每个碱基对应4个峰值。
                    # output_matrix是所有基因型在inference_window范围内的峰值数据
                    output_matrix[edit_proposal_idx][seq_index * 4 + color_index] = ep.trace_data[base_color_index]

                # normalize values
                # 标准化：每个碱基上的G、A、C、T的占比
                sum = np.sum(output_matrix[edit_proposal_idx][seq_index * 4:seq_index * 4 + 4])
                for color_index in range(4):
                    normalized_amt = output_matrix[edit_proposal_idx][seq_index * 4 + color_index] / sum
                    if np.isnan(normalized_amt * 100):
                        # 处理NaN（为什么会出现这个）
                        warnings.warn("空值")
                        output_matrix[edit_proposal_idx][seq_index * 4 + color_index] = 0
                    else:
                        output_matrix[edit_proposal_idx][seq_index * 4 + color_index] = normalized_amt * 100
        if self.verbose:
            print("Shape of coefficient matrix:", output_matrix.shape)
        # print(output_matrix)
        # 矩阵转置
        self.coefficient_matrix = output_matrix.T

    def _generate_outcomes_vector(self):
        # generate_b (outcomes) vector

        output_list = []
        edited_peak_values = self.edited_sample.get_peak_values()
        # print(edited_peak_values)
        inference_length = self.inference_window[1] - self.inference_window[0]
        # 诡异,是不是写反了   self.inference_window[0]:self.inference_window[1]
        self.inference_sequence = self.control_sample.primary_base_calls[
                                  self.inference_window[0]:self.inference_window[1]]

        last_good_ref_idx = -1
        for aln_tuple in self.alignment.alignment_pairs:
            ref_idx = aln_tuple[0]
            sample_idx = aln_tuple[1]
            if ref_idx is None:
                ref_idx = last_good_ref_idx
            else:
                last_good_ref_idx = ref_idx
            # ref_idx大于inference_window左端且output_list长度不能大于窗口大小 * 4,也就是填满之后就不给填了
            # 即在alignment_pairs中找到ref_idx在inference_window内的 下标对（控制组下标和实验组下标）
            if self.inference_window[0] <= ref_idx and len(output_list) < inference_length * 4:
                # 获取峰值
                for base in self.edited_sample.base_order:
                    # print(sample_idx)
                    val = edited_peak_values[base][sample_idx]
                    output_list.append(val)
        if self.verbose:
            print('------------------------------------------------------')
            print('Inference Sequence length: %d' % inference_length)
            print('Inference sequence', self.inference_sequence)
            print('Output vector : 4x%d (%d)' % (len(output_list) / 4, len(output_list)))

        self.output_vec = np.array(output_list)

    def normalize_observed_trace(self):
        """
        This function normalizes the observed Sanger sequencing trace for the sample.
        For each base position, it normalizes the signal for A, T, C, G to the total signal
        :return: vector with normalized numbers for the observed trace, b
        """
        b = self.output_vec

        b_prime = b.reshape(int(len(b) / 4), 4)

        row_sums = b_prime.sum(axis=1)
        normed = b_prime / row_sums[:, np.newaxis]

        b_normed = normed.reshape(int(len(b)))

        # if there any nan's due to dividing by zero, convert to zeros
        np.nan_to_num(b_normed, copy=False)
        return b_normed

    def simple_discordance_algorithm(self):
        """
        This is a simple algorithm that measures the discordance of the observed signal with the control signal (after the cut site).
        The proportion of edited sequences is then taken to be
        1-(signal conforming to control)
        take the minimum signal across the entire seq that corresponds to the control.
        the estimate of the edited proportion, could be an underestimate actually due to coincident bases
        :param norm_b:
        :return:

        """
        # discordances = self.discordances
        edited_discordances = self.results.edited_discordances
        control_discordances = self.results.control_discordances
        unexplained_discord_signal = []
        # print(ra['edited_discord'])

        for np_idx, edited_discord_signal in np.ndenumerate(edited_discordances):
            idx = np_idx[0]
            if self.guide_targets[0].cutsite < idx < self.inference_window[1]:
                control_discord_signal = control_discordances[idx]

                unexplained = edited_discord_signal - control_discord_signal

                if unexplained < 0:
                    unexplained = 0

                unexplained_discord_signal.append(unexplained)

        self.results.ice_d = min(100.,
                                 SangerAnalysis.iced_correction_factor * sum(unexplained_discord_signal) / float(
                                     len(unexplained_discord_signal)) * 100)

        self.results.max_unexp_discord = max(unexplained_discord_signal) * 100

    def infer_abundances(self, norm_b=False):
        """
        非负最小二乘法求每个基因型的丰度
        This uses nnls to solve for abundances of each edit proposal
        :param norm_b:
        :return:
        """
        # 基因型峰值（标准化）
        A = self.coefficient_matrix
        # 实验组峰值（未标准化）
        b = self.output_vec
        # import pdb; pdb.set_trace()
        if self.verbose:
            print("")
            print('NNLS input shapes')
            print("--------------------------------------------------------")
            print('A', A.shape)
            print('b', b.shape)

        # 默认为True，实验组数据标准化
        if norm_b:
            b_normed = self.normalize_observed_trace()
            b = b_normed
        # Solves argmin_x || Ax - b ||_2 for x>=0
        # A: columns (number of different indel possibilities)
        #  : rows (base calls (4*inference_length) + 1)
        # x: abundances of possibilities, has shape (A.cols, 1)
        # b: actual observed base calls, has shape (A.rows, 1)
        # nnls solves for x, or the abundances of each possible indel
        # xvals contains the inferred sequence abundances
        # rnorm is the residual || Ax-b ||_2

        try:

            # xvals, rnorm = nnls(A, b)
            # print(list(xvals))
            # compute the predicted signal
            # predicted = np.dot(A, xvals)

            # optional L1
            # 创建lasso模型
            lasso_model = linear_model.Lasso(alpha=0.8, positive=True)
            # print(lasso_model)
            lasso_model.fit(A, b)

            # 每一列都有一个系数（n行1列）n=列数（很多都是0，lasso模型的特点）
            xvals = lasso_model.coef_
            # print(list(xvals))
            # 列代表基因型，行代表GATC的峰值
            # 预测值，矩阵A点乘系数矩阵
            # 最后的结果是n行1列  n=行数
            # 相当于每个Base的预测值等于 每个基因型（每一列）在这个Base的峰值乘以这一列的系数 之和
            predicted = np.dot(A, xvals)
            # print(list(predicted))
            # print(list(b))
            # import csv
            # with open('C:/Users/41518/Desktop/data.csv','w',newline='') as aa:
            #     writer = csv.writer(aa)
            #     for i,j in zip(list(predicted),list(b)):
            #         writer.writerow([i,j])

        except Exception as e:

            raise type(e)(str(e) + ' A: ' + str(A.shape) + ' B: ' + str(b.shape))

        # calculate pearson's R
        (fit_r, p_val_2_tailed) = pearsonr(predicted, b)
        # fit_r 皮尔逊相关系数
        self.results.r_squared = fit_r ** 2
        print("R_SQUARED {}".format(self.results.r_squared))

        # 每种基因型的回归系数之和
        xtotal = xvals.sum()
        # here we normalize the relative abundance of each possible indel
        # 每种基因型的相对丰度
        for n, x_val in enumerate(xvals):
            self.proposals[n].x_abs = x_val
            self.proposals[n].x_rel = x_val / (1.0 * xtotal)

            if self.r_squared_correction:  # 默认值为True
                ##addition of (1-r_squared, or missing variance) to the no edit case
                if n == 0:
                    self.proposals[0].x_rel = x_val / (1.0 * xtotal) * self.results.r_squared
                ##edited cases
                else:
                    self.proposals[n].x_rel = x_val / (1.0 * xtotal) * self.results.r_squared
            # print(self.proposals[n].x_rel,self.proposals[n].bases_changed)
        new_array = round_percent([x.x_rel for x in self.proposals], self.results.r_squared)
        # np.sum(new_array)==r_squared
        # new_array=[np.floor(x.x_rel*100)/100 for x in self.proposals]

        for n, val in enumerate(new_array):
            self.proposals[n].x_rel = val
        self.results.r_squared = np.round(self.results.r_squared, 2)

    def analyze_and_rank(self):

        # 获取每个基因型的丰度和R2
        self.infer_abundances(norm_b=True)

        # 将基因型按照丰度大小排序
        sorted_by_contribution = sorted(self.proposals, key=lambda x: x.x_rel, reverse=True)
        # for i in sorted_by_contribution:
        #     print(i.x_rel)
        # create totals:

        aggregated_indel = defaultdict(float)
        hdr_percentage = 0
        ko_score = 0

        for entry in sorted_by_contribution:
            # 用于画分布图
            if entry.bases_changed not in aggregated_indel:
                # print(entry.bases_changed)
                # bases_changed是改变的碱基数量，如果不在aggregated_indel字典里，则添加
                aggregated_indel[entry.bases_changed] = 0
            aggregated_indel[entry.bases_changed] += entry.x_rel
            if entry.bases_changed == ProposalBase.HDR:
                # 这是啥？
                hdr_percentage += 100 * entry.x_rel

            # 如果bases_changed不是3的倍数或大于20，计算ko_score
            if entry.bases_changed != ProposalBase.HDR and (entry.bases_changed % 3 != 0 or abs(entry.bases_changed) > 20):
                ko_score += 100 * entry.x_rel
        # print(aggregated_indel)
        # print(aggregated_indel[0],self.results.r_squared)
        unedited_percent = aggregated_indel[0] + (1 - self.results.r_squared)

        editing_efficiency = (1.0 - unedited_percent)
        # print(editing_efficiency)
        if self.debug:
            print(aggregated_indel)
            print(sorted_by_contribution)
        self.results.ko_score = ko_score
        self.results.aggregate = aggregated_indel
        self.results.contribs = sorted_by_contribution
        self.results.edit_eff = editing_efficiency
        self.results.hdr_percentage = hdr_percentage

    #############################

    def analyze_sample(self):
        # 质量控制,找到平均质量分数大于40分的区域

        self.quality_check()
        # 找靶点
        self.find_targets()
        # print(self.guide_targets[0].cutsite)

        #
        if self.donor_odn:
            self.check_recombination()

        # 找比对窗口
        self.find_alignment_window()

        # 创建PairAlignment类，传入参数两个文件的序列
        alignment = PairAlignment(self.control_sample.primary_base_calls,
                                  self.edited_sample.primary_base_calls
                                  )

        self.alignment = alignment
        # 全局比对(调用align_all方法，改变自身属性)
        alignment.align_all()  # all_aligned_seqs和all_aligned_clustal属性改变
        # 局部比对（比对窗口）

        alignment.align_with_window(self.alignment_window)

        if alignment.has_alignment:
            pass
        else:
            raise Exception("No alignment found between control and edited sample")

        # 输出比对结果
        # aln_flag, msg = self.align_edited_to_control_sequence()
        aln_file = self.base_outputname + "all.txt"
        alignment.write_aln(alignment.all_aligned_clustal, to_file=aln_file)
        aln_file = self.base_outputname + "windowed.txt"
        alignment.write_aln(alignment.aln_clustal, to_file=aln_file)

        # write to jsons
        aln_json_file = self.base_outputname + "all.json"
        alignment.write_json(alignment.all_aligned_seqs, aln_json_file)
        aln_json_file = self.base_outputname + "windowed.json"
        alignment.write_json(alignment.aln_seqs, aln_json_file)

        # 找出所有可能的基因型
        self._generate_edit_proposals()

        # print(len(self.proposals))
        # 暂时忽略
        if self.donor_odn:
            aln_file = self.base_outputname + "donor.txt"
            self.donor_alignment.write_aln(self.donor_alignment.all_aligned_clustal, to_file=aln_file)
            aln_json_file = self.base_outputname + "donor.json"
            self.donor_alignment.write_json(self.donor_alignment.all_aligned_seqs, aln_json_file)

        print("analyzing {} number of edit proposals".format(len(self.proposals)))

        # 找出inference_window窗口
        self._calculate_inference_window()
        # 获取所有基因型的峰值数据（inference_window内）
        self._generate_coefficient_matrix()
        # 获取实验组的峰值数据（inference_window内）
        self._generate_outcomes_vector()

        self.analyze_and_rank()

        self.calculate_discordance()
        self.simple_discordance_algorithm()

        # output
        indel_file = self.base_outputname + "indel.json"
        trace_file = self.base_outputname + "trace.json"

        generate_discordance_indel_files(self, self.results, to_file=indel_file)
        generate_trace_files(self, to_file=trace_file)

        write_individual_contribs(self, to_file=self.base_outputname + "contribs.txt")
        write_contribs_json(self, self.base_outputname + "contribs.json")

        if self.allprops:
            write_all_proposals_json(self, self.base_outputname + "allproposals.json")
