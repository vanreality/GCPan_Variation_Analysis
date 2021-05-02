import argparse
import os


def get_query(record):
    return record.strip().split("\t")[0]


def get_search(record):
    return record.strip().split("\t")[1]


def main(args):
    with open("185samples", 'r') as f:
        samples = f.readlines()

    root = "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result"
    manta_dir = root + "/manta"
    blastn_dir = root + "/blastn"
    mapped_novel = {}
    mapped_ref = {}
    mapped_none = {}
    for sample in samples:
        sample = sample.strip()
        for sv_type in ["germline", "somatic"]:
            path = blastn_dir + "/" + sample + "/" + sv_type + "SV.blastn.out.filtered"
            if os.path.exists(path):
                with open(path) as f:
                    for line in f.readlines():
                        if line[0] == "#":
                            continue
                        if "GC" in get_search(line):
                            mapped_novel[sample + ":" + get_query(line)] = ""
                        else:
                            mapped_ref[sample + ":" + get_query(line)] = ""
            path = manta_dir + "/" + sample + "/result/results/variants/" + sv_type + "SV.vcf.fasta"
            if os.path.exists(path):
                with open(path) as f:
                    ins_id = ""
                    for line in f.readlines():
                        if line[0] == ">":
                            ins_id = sample + ":" + line[1:].split("\t")[0].split(" ")[0]
                        else:
                            if ins_id in mapped_novel:
                                mapped_novel[ins_id] = ">" + ins_id.split(".")[0] + "\n" + line
                            elif ins_id in mapped_ref:
                                mapped_ref[ins_id] = ">" + ins_id.split(".")[0] + "\n" + line
                            else:
                                mapped_none[ins_id] = ">" + ins_id.split(".")[0] + "\n" + line

    with open(args.o + "mapped_to_novel.fasta", 'w') as f:
        for ins_id in mapped_novel:
            f.write(mapped_novel[ins_id])

    with open(args.o + "mapped_to_ref.fasta", 'w') as f:
        for ins_id in mapped_ref:
            f.write(mapped_ref[ins_id])

    with open(args.o + "mapped_to_none.fasta", 'w') as f:
        for ins_id in mapped_none:
            f.write(mapped_none[ins_id])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", help="output file prefix")

    main(parser.parse_args())
