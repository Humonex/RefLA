class Scaffold:
    def __init__(self):
        self.contigs = []
        self.overlaps = {}

    def add_contig(self, contig_id, sequence):
        self.contigs.append((contig_id, sequence))

    def add_overlap(self, contig1_id, contig2_id, overlap_length):
        self.overlaps[(contig1_id, contig2_id)] = overlap_length

    def scaffold_contigs(self):
        # Sort contigs based on the number of overlaps
        sorted_contigs = sorted(self.contigs, key=lambda x: len(self._get_overlaps(x[0])), reverse=True)

        # Start scaffolding from contig with the highest number of overlaps
        scaffold_sequence = sorted_contigs[0][1]
        used_contigs = {sorted_contigs[0][0]}

        # Loop through the sorted contigs and join them based on overlaps
        for contig_id, sequence in sorted_contigs[1:]:
            if contig_id not in used_contigs:
                for u_contig_id in used_contigs:
                    overlap = self.overlaps.get((u_contig_id, contig_id)) or self.overlaps.get((contig_id, u_contig_id))
                    if overlap:
                        # Join the sequences based on overlap
                        index = scaffold_sequence.find(sequence[:overlap])
                        if index != -1:
                            scaffold_sequence = scaffold_sequence[:index] + sequence
                            used_contigs.add(contig_id)
                            break

        return scaffold_sequence

    def _get_overlaps(self, contig_id):
        return [k for k in self.overlaps if contig_id in k]

# Example usage
scaffolder = Scaffold()
scaffolder.add_contig('CM008963.1:3622156-3766256', 'sequence1')
scaffolder.add_contig('CM008963.1:2863659-2991359', 'sequence2')
scaffolder.add_overlap('CM008963.1:3622156-3766256', 'CM008963.1:2863659-2991359', 50)

# Add more contigs and overlaps as needed...

scaffolded_sequence = scaffolder.scaffold_contigs()
print(scaffolded_sequence)

from Bio.Seq import Seq
from collections import defaultdict
 
seq_dict = {
    'seq1': Seq('AACGTTGGAA'),
    'seq2': Seq('TTGCCGTTAA'),
    'seq3': Seq('GGCTTAAACG'),
    'seq4': Seq('TTAAGGCCCG')
}
 
bridge_sequences = ['TTGG', 'GCCG', 'GGCT']
 
def assign_id(sequence_dict, bridge_seq):
    id_bridge_dict = defaultdict(list)
 
    for key, sequence in sequence_dict.items():
        for i in range(len(sequence) - len(bridge_seq[0]) + 1):
            sub_seq = str(sequence[i:i+len(bridge_seq[0])])
            if sub_seq in bridge_seq:
                id_bridge_dict[key].append(sub_seq)
    return dict(id_bridge_dict)
 
result = assign_id(seq_dict, bridge_sequences)
print(result)