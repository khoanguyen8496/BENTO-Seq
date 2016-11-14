import pandas as pd
import json

def parser(fevents, bio_name_hash="bio_symbol_ensg.txt", fout="base_data.json"):
    """
    Prepare JSON file for all splicing events
    :param fevents:
    :return: JSON object of all events
    """

    hash = json.load(open(bio_name_hash))
    df = pd.read_csv(fevents, header=0, sep="\t")
    uKeys = list(df["id"])
    genes_names = [x.split(":")[2] for x in uKeys]
    triplet_ids = map(int, [x.split(":")[4] for x in uKeys])
    gene_ensgs = map(lambda x: hash.get(x, x), genes_names)

    c1_s = map(lambda x: map(int, x.split(":"))[0], df["Up-stream Exon"])
    c1_e = map(lambda x: map(int, x.split(":"))[1], df["Up-stream Exon"])
    c2_s = map(lambda x: map(int, x.split(":"))[0], df["Down-stream Exon"])
    c2_e = map(lambda x: map(int, x.split(":"))[1], df["Down-stream Exon"])
    a_s = map(lambda x: map(int, x.split(":"))[0], df["Cassette Exon"])
    a_e = map(lambda x: map(int, x.split(":"))[1], df["Cassette Exon"])

    df["ensg"] = gene_ensgs
    df["triplet"] = triplet_ids
    df["c1_e"] = c1_e
    df["c2_e"] = c2_e
    df["c2_e"] = c2_e
    df["c1_s"] = c1_s
    df["a_S"] = a_s
    df["a_s"] = a_s
    df["a_e"] = a_e
    df["c2_s"] = c2_s

    base_data = {x:{"gene": df[x]["ensg"],
                "strand": df[x]["strand"],
                "triplet": df[x]["triplet"],
                "C1_start":df[x]["c1_s"],
                "C1_end": df[x]["c1_e"],
                "A_start": df[x]["a_s"],
                "A_end": df[x]["a_e"],
                "C2_start": df[x]["c2_s"],
                "C2_end": df[x]["c2_e"],
                } for x in df.keys()}

    json.dump(base_data, open("base_data.json", "w"))
    pass