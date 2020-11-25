import sys


class Validation:
    def __init__(self, reference_file_name, prediction_file_name):
        self._reference_file_name = reference_file_name
        self._prediction_file_name = prediction_file_name
        self._ref_pos_to_mutation = {}
        self._pred_pos_to_mutation = {}
        self.fill_pos_to_mutation()  # Fill reference
        self.fill_pos_to_mutation(False)  # Fill predictions

    def fill_pos_to_mutation(self, ref=True):
        if ref:
            pos_to_mutation = self._ref_pos_to_mutation
            file_name = self._reference_file_name
        else:
            pos_to_mutation = self._pred_pos_to_mutation
            file_name = self._prediction_file_name

        with open(file_name, "r") as file:
            for line in file:
                if line[0] == "#":
                    continue
                # line example:
                # 126     A       T       2
                pos = int(line.split()[0])
                alt = line.split()[2]
                if pos not in pos_to_mutation:
                    pos_to_mutation[pos] = set()
                pos_to_mutation[pos].add(alt)

    def print_precision_recall(self):
        nb_variants = 0
        nb_TP_variants = 0
        nb_FN_variants = 0
        for pos, alts in self._ref_pos_to_mutation.items():
            nb_variants += len(alts)
            if pos in self._pred_pos_to_mutation:
                for alt in alts:
                    if alt in self._pred_pos_to_mutation[pos]:
                        nb_TP_variants += 1
                    else:
                        nb_FN_variants += 1
            else:
                nb_FN_variants += len(alts)

        nb_FP_variants = 0
        for pos, alts in self._pred_pos_to_mutation.items():
            if pos in self._ref_pos_to_mutation:
                for alt in alts:
                    if alt not in self._ref_pos_to_mutation[pos]:
                        nb_FP_variants += 1
            else:
                nb_FP_variants += len(alts)

        print(f"Nb Variants {nb_variants} in truth")
        print(f"Nb TP {nb_TP_variants} in truth & in pred")
        print(f"Nb FP {nb_FP_variants} not in truth & in pred")
        print(f"Nb FN {nb_FN_variants} in truth & not in pred")
        print(f"Precision {round(100*nb_TP_variants/float(nb_TP_variants+nb_FP_variants),1)} %")
        print(f"Recall {round(100*nb_TP_variants/float(nb_TP_variants+nb_FN_variants),1)} %")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} truth.vcf to_test.vcf")
        exit()
    validation = Validation(sys.argv[1], sys.argv[2])
    validation.print_precision_recall()
