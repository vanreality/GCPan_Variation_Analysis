import argparse
import os


def cal_cov_len(data):
    sorted_by_lower_bound = sorted(data, key=lambda tup: tup[0])
    union = []
    cov_len = 0

    for higher in sorted_by_lower_bound:
        if not union:
            union.append(higher)
        else:
            lower = union[-1]
            if higher[0] > lower[1]:
                union.append(higher)
            else:
                union[-1] = (lower[0], higher[1])

    for interval in union:
        cov_len += interval[1] - interval[0] + 1

    return cov_len


def main(args):
    s = open(args.s, 'r')
    l = open(args.l, 'r')
    n = open(args.n, 'r')
    r = open(args.i + "/../coverage", 'w')
    # r = open(args.i + "/../coverage.unfiltered", 'w')  # unfiltered coverage

    cds_len = {}
    result = {}
    cov_len = {}
    index = 0

    for line in l.readlines():
        line = line.split("\t")
        cds_len[line[0]] = line[2].strip()

    for sample in s.readlines():
        # if os.path.exists(args.i + "/" + sample.strip()):
        #     with open(args.i + "/" + sample.strip()) as s:
        if os.path.exists(args.i + "/" + sample.strip()):
            with open(args.i + "/" + sample.strip()) as s:
                for line in s.readlines():
                    line = line.split("\t")
                    gene = line[7]

                    start = max(int(line[1]), int(line[5]))
                    end = min(int(line[2]), int(line[6]))

                    if gene not in cov_len:
                        cov_len[gene] = []
                    else:
                        cov_len[gene].append([start, end])

            for gene in cov_len:
                cov = 1 - cal_cov_len(cov_len[gene]) / int(cds_len[gene])

                if gene not in result:
                    result[gene] = [1] * 185

                result[gene][index] = round(cov, 4)

            cov_len = {}
            index += 1

    sample_names = "Gene"
    for sample in n.readlines():
        sample = sample.split("\t")
        sample_names += "\t" + sample[1].strip()
    r.write(sample_names + "\n")

    for gene in result:
        line = gene
        for cov in result[gene]:
            line += "\t" + str(cov)
        r.write(line + "\n")

    s.close()
    r.close()
    n.close()
    r.close()


if __name__ == "__main__":
    root = "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="intersect bed file directory",
                        default=root + "result/compareSVPAV/svaba.germline.sv")
    parser.add_argument("-s", help="input samples",
                        default=root + "script/analyzeSV/185samples")
    parser.add_argument("-n", help="sample name transfer table",
                        default=root + "data/PAV/normal.individual.info")
    parser.add_argument("-l", help="cds length table",
                        default=root + "data/PAV/cds.length")
    args = parser.parse_args()

    main(args)
