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
    if get_score(a) > get_score(b):
        return a
    elif get_score(a) == get_score(b):
        if get_indent(a) > get_indent(b):
            return a
        elif get_indent(a) == get_indent(b):
            if get_length(a) > get_length(b):
                return a
            else:
                return b
        else:
            return b
    else:
        return b


def main(arg):
    best_hit = {}
    with open(arg.i, 'r') as i:
        for line in i.readlines():
            query = get_query(line)
            if query not in best_hit:
                if "GC" in get_search(line):
                    best_hit[query] = []
                    best_hit[query].append("")
                    best_hit[query].append(line)
                else:
                    best_hit[query] = []
                    best_hit[query].append(line)
                    best_hit[query].append("")
            else:
                if "GC" in get_search(line):
                    if best_hit[query][1] == "":
                        best_hit[query][1] = line
                    else:
                        best_hit[query][1] = compare(best_hit[query][1], line)
                else:
                    if best_hit[query][0] == "":
                        best_hit[query][0] = line
                    else:
                        best_hit[query][0] = compare(best_hit[query][0], line)

    with open(arg.i + ".filtered", 'w') as o:
        for query in best_hit:
            if best_hit[query][1] != "":
                o.write(best_hit[query][0])
                o.write(best_hit[query][1])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input blastn result file name")
    args = parser.parse_args()

    main(args)
