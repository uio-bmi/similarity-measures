import pandas
from scipy.stats import rankdata
from scipy.stats.stats import pearsonr
import numpy as np


def create_dict(filename):
    scoretable = pandas.read_csv(filename, sep="\t", header=2, index_col=0)
    ijs = [(i, j) for i in range(scoretable.shape[0])
           for j in range(i+1, scoretable.shape[1])]
    keys = [scoretable.index[i] + scoretable.columns[j] for i, j in ijs]
    values = [scoretable.values[i, j] for i, j in ijs]
    ranks = rankdata(values, method="ordinal")
    # print(filename, max(ranks))
    mv = max(values)
    # print(mv, sum(1 for v in values if v == mv))
    return dict(zip(keys, ranks))


def get_correlation(keys, dicts, metrics):
    values = {metric: [d[key] for key in keys]
              for metric, d in zip(metrics, dicts)}
    pairs = [(metric_a, metric_b)
             for i, metric_a in enumerate(metrics)
             for metric_b in metrics[i+1:]]
    for a, b in pairs:
        print(a, b, pearsonr(values[a], values[b]))


def get_top_overlap(keys, dicts, metrics):
    values = {metric: np.array([d[key] for key in keys])
              for metric, d in zip(metrics, dicts)}
    # print({k: (np.max(v), len(v)) for k, v in values.items()})
    top_values = {metric: np.argsort(v)[-100:]
                  for metric, v in values.items()}
    top_keys = {metric: {keys[i] for i in tops}
                for metric, tops in top_values.items()}
    pairs = [(metric_a, metric_b)
             for i, metric_a in enumerate(metrics)
             for metric_b in metrics[i+1:]]
    # print({key: len(v) for key, v in top_keys.items()})
    for a, b in pairs:
        print(a, b, len(top_keys[a] & top_keys[b]))

    print("All", len(top_keys["forbes"] & top_keys["jaccard"] & top_keys["tetra"]))


def main():
    metrics = ["jaccard", "forbes", "tetra"]
    dicts = [create_dict(metric+"/scores.txt") for metric in metrics]
    keys = list(dicts[0].keys())
    print("Correleation")
    get_correlation(keys, dicts, metrics)
    print("Top 100 Overlap")
    get_top_overlap(keys, dicts, metrics)
    rows = [[key]+[str(d[key]) for d in dicts] for key in keys]
    outfile = open("all_ranks.csv", "w")
    outfile.write("Pair\t%s\n" % "\t".join(metrics))
    for row in rows:
        outfile.write("\t".join(row)+"\n")

if __name__ == "__main__":
    main()
