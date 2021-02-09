import re
class Donor:
    def __init__(self,donor):
        self.homologous_arm = 20
        self.donor_seq = donor

        self.aligned_seq = None
        self.aligned_pairs_index = []

        self.genotypes = None
        self.change_bases = None

        self.donor_pos = None


    def make_donor_index(self):
        # make index
        control_index = []
        e_index = 0
        donor_index = []
        d_index = 0
        first_align_index = re.search('[A-Z]', self.aligned_seq[1]).start()
        last_align_index = re.search('[A-Z\-]+[A-Z]', self.aligned_seq[1]).end()
        self.donor_pos = (first_align_index, last_align_index)
        self.aligned_seq[1] = self.aligned_seq[0][:first_align_index] + self.aligned_seq[1][first_align_index:last_align_index] + self.aligned_seq[0][last_align_index:]
        for i in range(len(self.aligned_seq[0])):
            if self.aligned_seq[0][i] != '-':
                control_index.append(e_index)
                e_index += 1
            else:
                control_index.append(None)

            if self.aligned_seq[1][i] != '-':
                donor_index.append(d_index)
                d_index += 1
            else:
                donor_index.append(None)


        for id1, id2 in zip(control_index, donor_index):
            self.aligned_pairs_index.append([id1, id2])
