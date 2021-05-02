import pandas as pd


root = "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result"
blastn_dir = root + "/blastn"
trf_dir = root + "/trf"
threshold = 0.5

with open("185samples", 'r') as f:
    samples = {}
    for sample in f.readlines():
        samples[sample.split(".")[0]] = sample.strip()

res = pd.read_csv(blastn_dir + "/summary.0.9.csv", header=1, index_col=0)
for mapped_type in ["novel", "ref", "none"]:
    with open(trf_dir + "/mapped_to_" + mapped_type + ".fasta.2.7.7.80.10.50.500.mask", 'r') as f:
        masked = {}
        seq = ""
        sample = ""
        for line in f.readlines():
            if line[0] == ">":
                if seq != "":
                    if (seq.count("N") / len(seq)) >= threshold:
                        masked[samples[sample]] += 1
                seq = ""
                sample = line[1:].strip()
                if samples[sample] not in masked:
                    masked[samples[sample]] = 0
            else:
                seq += line.strip()

    for sample in masked:
        res.loc["mapped_" + mapped_type, sample] -= masked[sample]

res.to_csv(blastn_dir + "/summary.masked.csv")
