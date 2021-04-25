import argparse
import os
import pandas as pd


def union_covered(data):
    sorted_by_lower_bound = sorted(data, key=lambda tup: tup[0])
    union = []

    for higher in sorted_by_lower_bound:
        if not union:
            union.append(higher)
        else:
            lower = union[-1]
            if higher[0] > lower[1]:
                union.append(higher)
            else:
                union[-1] = (lower[0], max(lower[1], higher[1]))

    return union


def count_covered_length(data):
    length = 0
    for i in data:
        length += i[1] - i[0]
    return length


def main(args):
    root = "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/blastn"

    with open(args.i, 'r') as s:
        samples = s.readlines()

    # summarize INS aligned to neither reference nor novel sequence, and INS aligned to novel sequence
    # and calculate total length in novel sequence covered by SV INS
    cols = []
    novel_sequence_covered = {}
    for sample in samples:
        cols.append(sample.strip())
    summary = pd.DataFrame(index=["mapped_novel", "mapped_ref", "no_mapped", "total_ins"], columns=cols)
    for sample in samples:
        sample = sample.strip()
        mapped_novel, mapped_ref, no_mapped, total_ins = [0] * 4
        for var_type in ["germline", "somatic"]:
            path = root + "/" + sample + "/" + var_type + "SV.blastn.out.filtered"
            if os.path.exists(path):
                with open(path, 'r') as f:
                    lines = f.readlines()
                    total_ins += int(lines[0].strip().split(":")[1])
                    mapped_novel += int(lines[1].strip().split(":")[1])
                    mapped_ref += int(lines[2].strip().split(":")[1])
                    for line in lines[3:]:
                        line = line.strip().split("\t")
                        if line[1] not in novel_sequence_covered:
                            novel_sequence_covered[line[1]] = []
                        lower_bound = min(int(line[-5]), int(line[-4]))
                        higher_bound = max(int(line[-5]), int(line[-4]))
                        novel_sequence_covered[line[1]].append([lower_bound, higher_bound])
        no_mapped = total_ins - mapped_ref - mapped_novel
        summary[sample] = [mapped_novel, mapped_ref, no_mapped, total_ins]

    covered_length = 0
    for seq in novel_sequence_covered:
        novel_sequence_covered[seq] = union_covered(novel_sequence_covered[seq])
        covered_length += count_covered_length(novel_sequence_covered[seq])

    # output summary
    summary.to_csv(args.o)
    with open(args.o, 'r+') as o:
        content = o.read()
        o.seek(0, 0)
        o.write("# Total length covered by SV INS : " + str(covered_length) + "\n" + content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input samples names")
    parser.add_argument("-o", help="output file name")

    main(parser.parse_args())
