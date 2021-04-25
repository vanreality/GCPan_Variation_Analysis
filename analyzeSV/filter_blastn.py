import argparse


def get_query(record):
    return record.strip().split("\t")[0]


def get_search(record):
    return record.strip().split("\t")[1]


def get_indent(record):
    return float('%.1f' % float(record.strip().split("\t")[2]))


def get_score(record):
    return float('%.1f' % float(record.strip().split("\t")[11]))


def get_length(record):
    return int(record.strip().split("\t")[3])


def compare(a, b):
    if get_indent(a) * get_length(a) > get_indent(b) * get_length(b):
        return a
    elif get_indent(a) * get_length(a) == get_indent(b) * get_length(b):
        if get_score(a) > get_score(b):
            return a
        else:
            return b
    else:
        return b
    # if get_score(a) > get_score(b):
    #     return a
    # elif get_score(a) == get_score(b):
    #     if get_indent(a) > get_indent(b):
    #         return a
    #     elif get_indent(a) == get_indent(b):
    #         if get_length(a) > get_length(b):
    #             return a
    #         else:
    #             return b
    #     else:
    #         return b
    # else:
    #     return b


def main(arg):
    ins_length = {}
    with open(args.f, 'r') as f:
        for line in f.readlines():
            if line[0] == ">":
                ins_id = line[1:].split("\t")[0].split(" ")[0]
                ins_length[ins_id] = 0
            else:
                ins_length[ins_id] = len(line.strip())

    best_hit = {}
    with open(arg.i, 'r') as i:
        for line in i.readlines():
            query = get_query(line)
            if get_length(line) / ins_length[query] >= 0.9:
                if query not in best_hit:
                    best_hit[query] = line
                else:
                    best_hit[query] = compare(best_hit[query], line)

    if best_hit != {}:
        with open(arg.i + ".filtered", 'w') as o:
            # filter at least 90% length coverage
            map_to_novel = 0
            map_to_ref = 0
            for query in best_hit:
                if "GC" in get_search(best_hit[query]):
                    map_to_novel += 1
                else:
                    map_to_ref += 1
            o.write("# Total INS sequence number:" + str(len(ins_length)) + "\n")
            o.write("# INS mapped to novel sequence:" + str(map_to_novel) + "\n")
            o.write("# INS mapped to reference:" + str(map_to_ref) + "\n")
            for query in best_hit:
                o.write(best_hit[query].strip() + "\t" + str(ins_length[query]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input blastn result file name")
    parser.add_argument("-f", help="input manta fasta file name")
    args = parser.parse_args()

    main(args)
