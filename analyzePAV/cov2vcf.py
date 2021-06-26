import argparse
import pandas as pd


def main(args):
    with open(args.v, 'rt') as h:
        vcf = h.readlines()
    bed = pd.read_table(args.b, header=None)
    with open(args.c, 'rt') as c:
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        for i in c.readline().strip().split("\t")[1:]:
            header += "\t" + i
        vcf.append(header + "\n")

        for line in c.readlines():
            line = line.strip().split("\t")
            tmp = bed.loc[bed[3] == line[0]].values[0]
            record = ""
            record += (tmp[0] + "\t" + str(tmp[1]) +
                       "\t" + tmp[3] + "\t<PRE>\t<ABS>\t99\tPASS\t" + str(tmp[5]) + "\tGT")
            for i in line[1:]:
                record += "\t0/0" if float(i) > args.t else "\t1/1"
            vcf.append(record + "\n")

    with open(args.o, 'w') as o:
        for line in vcf:
            o.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="input coverage file name")
    parser.add_argument("-v", help="input header name",
                        default="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/PAV/vcf_header.txt")
    parser.add_argument("-t", help="coverage threshold",
                        default=0.8)
    parser.add_argument("-b", help="gene bed file name",
                        default="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/PAV/pan.bed")
    parser.add_argument("-o", help="output vcf file name")

    main(parser.parse_args())
