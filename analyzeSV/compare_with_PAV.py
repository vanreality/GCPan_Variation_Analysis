import argparse


def main(args):
    s = open(args.s, "r")
    p = open(args.p, "r")
    o = open(args.s + "." + str(args.t) + ".count.csv", "w")
    r = open(args.s + "." + str(args.t) + ".csv", "w")
    g = open(args.s + "." + str(args.t) + ".gene.csv", "w")

    sv = {}
    output_count = {}
    output_compare = {}
    output_pav_specific_gene = {}

    for line in s.readlines():
        line = line.split("\t")
        if line[0] == "Gene":
            continue
        cov = []
        for i in range(1, len(line)):
            if float(line[i].strip()) < args.t:
                cov.append(0)
            else:
                cov.append(1)
        sv[line[0]] = cov

    sample_number = len(cov)
    title = ""

    for line in p.readlines():
        line = line.split("\t")
        if line[0] == "Gene":
            for tmp in line:
                title += tmp + ","
            title = title[0:-1]
            continue

        cov = []
        count = 0

        if line[0] not in sv:
            for i in range(sample_number):
                if float(line[i + 1].strip()) < args.t:
                    cov.append(3)
                    count += 1
                else:
                    cov.append(1)
            output_compare[line[0]] = cov
            if count >= 1:
                output_pav_specific_gene[line[0]] = cov
            continue

        for i in range(sample_number):
            if float(line[i + 1].strip()) < args.t:
                if sv[line[0]][i] == 0:
                    cov.append(0)  # SV absence
                else:
                    cov.append(3)  # noSV absence
                    count += 1
            else:
                if sv[line[0]][i] == 0:
                    cov.append(2)  # SV presence
                else:
                    cov.append(1)  # noSV presence

        output_compare[line[0]] = cov

        if count >= 1:
            output_pav_specific_gene[line[0]] = cov

    output_count[0] = [0] * sample_number
    output_count[1] = [0] * sample_number
    output_count[2] = [0] * sample_number
    output_count[3] = [0] * sample_number

    for i in output_compare:
        for j in range(sample_number):
            try:
                output_count[output_compare[i][j]][j] += 1
            except Exception as e:
                continue

    o.write(title)
    compare_type = ["SV_Absence", "NoSV_Presence", "SV_Presence", "NoSV_Absence"]
    for i in range(4):
        tmp = compare_type[i]
        for j in range(sample_number):
            tmp += "," + str(output_count[i][j])
        o.write(tmp + "\n")

    r.write(title)
    for gene in output_compare:
        tmp = gene
        for i in range(sample_number):
            tmp += "," + str(output_compare[gene][i])
        r.write(tmp + "\n")

    g.write(title)
    for gene in output_pav_specific_gene:
        tmp = gene
        for i in range(sample_number):
            tmp += "," + str(output_pav_specific_gene[gene][i])
        g.write(tmp + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help="input SV coverage file")
    parser.add_argument("-p", help="input PAV coverage file")
    parser.add_argument("-t", default=0.8, type=float, help="coverage threshold")
    args = parser.parse_args()

    main(args)
