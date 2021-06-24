import argparse


def main(args):
    with open(args.h, 'rt') as h:
        header = h.readlines()
    with open(args.c, 'rt') as c:
        format = c.readline()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="input coverage file name")
    parser.add_argument("-h", help="input header name", default="vcf_header.txt")
    parser.add_argument("-t", help="coverage threshold", default=0.8)
    parser.add_argument("-o", help="output vcf file name")

    main(parser.parse_args())
