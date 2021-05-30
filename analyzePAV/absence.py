import argparse


def main(args):
    print("Hello")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input file name")

    main(parser.parse_args())
