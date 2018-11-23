from scoretable import ScoreTable, TrackInfo, strip_title
from math import log


def generate_data(score_table, track_info):
    rows = []
    for title, name in zip(score_table._table.columns,
                           score_table.tf_types):
        score_column = score_table._table[title]
        for title2, name2 in zip(
                score_column.index, score_table.tf_types):
            if name < name2 or strip_title(title) == strip_title(title2):
                continue
            size1 = track_info.get_bp_covered(title)
            size2 = track_info.get_bp_covered(title2)
            row = [name + "_" + name2,
                   size1+size2, log(size1)+log(size2), score_column[title2]]
            rows.append(row)
    return rows


def generate_data2(score_table, track_info):
    rows = []
    for title, name in zip(score_table._table.columns,
                           score_table.tf_types):
        print("|" + strip_title(title) + "|")
        if strip_title(title).strip() == "MACS_ENCSR000DZM_GM12878_STAT1 MACS_ENCSR000DZM_GM12878_STAT1":
            print("Skipping MACS_ENCSR000DZM_GM12878_STAT1 MACS_ENCSR000DZM_GM12878_STAT1")
            continue
        score_column = score_table._table[title]
        for title2, name2 in zip(
                score_column.index, score_table.tf_types):
            if strip_title(title) == strip_title(title2):
                continue
            size1 = track_info.get_bp_covered(title)
            size2 = track_info.get_bp_covered(title2)
            row = [strip_title(title), name, size1,
                   strip_title(title2), name2, size2,
                   score_column[title2]]
            rows.append(row)
    return rows


def main(score_table_name, track_info_name):
    score_table = ScoreTable(score_table_name)
    track_info = TrackInfo(track_info_name)
    return generate_data2(score_table, track_info)

if __name__ == "__main__":
    import sys
    rows = main(sys.argv[1], sys.argv[2])
    text_rows = ([str(e) for e in row] for row in rows)
    with open(sys.argv[3], "w") as f:
        f.write("\t".join(
            ["title", "name", "size", "title2",
             "name2", "size2", "score"])+"\n")
        # f.write("\t".join(["group", "size", "logsize", "score"]) + "\n")
        f.writelines("\t".join(row)+"\n" for row in text_rows)
