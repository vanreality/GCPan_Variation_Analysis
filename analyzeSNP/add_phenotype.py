import argparse
import pandas as pd


def main(args):
    phenotype = pd.read_table(args.p)
    fam = pd.read_csv(args.f, sep=args.s, header=None)

    for index, row in phenotype.iterrows():
        fam.loc[fam[1] == row["Normal"], 4] = row["Gender"]
        fam.loc[fam[1] == row["Normal"], 5] = row[args.t]

    fam.to_csv(args.o, header=False, index=False, sep=args.s)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="input phenotype file")
    parser.add_argument("-f", help="input fam file")
    parser.add_argument("-t", help="selected phenotype")
    parser.add_argument("-s", help="sep", default="\t")
    parser.add_argument("-o", help="output file")

    main(parser.parse_args())
