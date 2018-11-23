import pandas
import statsmodels.api as sm
# from sklearn.preprocessing import LabelEncoder
from statsmodels.regression.linear_model import WLS
from sklearn import linear_model
import matplotlib.pyplot as plt
# from  matplotlib.axes import Axes
import numpy as np


def strip_title(name):
    return name.split("(")[0].strip()


def between_sd(group_codes, ys):
    ys = ys[:, None]
    ns = np.sum(group_codes, axis=0)
    sums = np.sum(group_codes*ys, axis=0)
    means = sums/ns
    grand_mean = np.sum(ys)/ys.size
    var = np.sum(ns*(means-grand_mean)**2)
    return np.sqrt(var/(ys.size-1))


class Analysis:
    def __init__(self, score_table, track_info):
        self._score_table = score_table
        self._track_info = track_info

    def plot(self):
        intercepts = np.sum(self.model.coef_[:-1]*self.xs[:, :-1], 1)
        t_normalized = self.ys  #  - intercepts
        normalized = self._inverse_y(t_normalized)
        plt.scatter(self.xs[:, -1], normalized)
        plt.xlabel(self.xlabel)
        plt.ylabel("Score")
        plt.title(self.title)
        np.savetxt(
            self.title+"_scatter.tsv", np.vstack((self.xs[:, -1], normalized)),
            delimiter="\t")
        plt.savefig(self.title+"_scoresvsize.png")
        plt.close()

        plt.scatter(self.xs[:, -1], t_normalized)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.savefig(self.title + "_t_scoresvsize.png")
        plt.close()
        np.savetxt(self.title + "_t_scatter.tsv", np.vstack((self.xs[:, -1], t_normalized)),
                   delimiter="\t")
        # end_points = np.linspace(np.min(self.xs[:, -1]),
        #                         np.max(self.xs[:, -1]), 50)
        #
        # line = self.model.intercept_ + self.model.coef_[-1]*end_points
        # line_p = line+self.sd
        # line_n = line-self.sd
        # line = self._inverse_y(line)
        # line_p = self._inverse_y(line_p)
        # line_n = self._inverse_y(line_n)
        # plt.plot(end_points, line, color="g")
        # np.savetxt("green_line.tsv", np.vstack((end_points, line)), delimiter="\t")
        # plt.plot(end_points, line_n, color="r")
        # np.savetxt("red_line.tsv", np.vstack((end_points, line_p)), delimiter="\t")
        # plt.plot(end_points, line_p, color="r")
        # np.savetxt("red_line2.tsv", np.vstack((end_points, line_p)), delimiter="\t")
        # plt.savefig("scoresvsize.png")
        # plt.close()

    def generate_data(self):
        self.pair_ids = []
        self.scores = []
        self.size_measures = []
        self.internal_idxs = []
        i = 0
        for title, name in zip(self._score_table._table.columns,
                               self._score_table.tf_types):
            score_column = self._score_table._table[title]
            for title2, name2 in zip(
                    score_column.index, self._score_table.tf_types):
                if name < name2 or strip_title(title) == strip_title(title2):
                    continue
                self.pair_ids.append(name+"_"+name2)
                self.size_measures.append(self.size_measure(
                    self._track_info.get_bp_covered(title),
                    self._track_info.get_bp_covered(title2)))
                self.scores.append(score_column[title2])
                if name == name2:
                    self.internal_idxs.append(i)
                i += 1
        self.encode()

    def encode(self):
        encoded_pairs = self.encode_data(self.pair_ids)
        self.xs = np.zeros((len(self.scores), encoded_pairs.shape[1]+1))
        self.xs[:, :-1] = encoded_pairs
        self.xs[:, -1] = self.size_measures
        # self.xs[self.internal_idxs, -1] = 1
        self.ys = self._transform_y(np.array(self.scores))
        self.sd = between_sd(self.xs[:, :-1], self.ys)

    def encode_data(self, data):
        values = list(set(data))
        lookup = {v: j for j, v in enumerate(values[1:])}
        enc_data = np.zeros((len(data), len(values)-1))
        for i, val in enumerate(data):
            if val not in lookup:
                continue
            enc_data[i, lookup[val]] = 1
        self.values = values
        self.lookup = lookup

        return enc_data

    def analyze(self):
        model = linear_model.LinearRegression()
        model.fit(self.xs, self.ys)
        self.model = model

    def _transform_y(self, ys):
        self.ylabel = "score"
        return ys

    def size_measure(self, size1, size2):
        self.xlabel = "log(size1) + log(size2)"
        return np.log(size1*size2+1)

    def _size_measure(self, size1, size2):
        self.xlabel = "size1+size2"
        assert size1 >= 0 and size2 >= 0
        return size1+size2


class JaccardAnalysis(Analysis):
    title = "Jaccard"

    def _size_measure(self, size1, size2):
        self.xlabel = "log(size1) + log(size2)"
        return np.log(size1*size2+1)

    def __transform_y(self, ys):
        self.ylabel = "logit(score)"
        assert np.all(ys >= 0) and np.all(ys <= 1)
        return ys  # np.log(ys+0.000001)  # (ys+0.0001)/(1.0001-ys))

    def _inverse_y(self, ys):
        return (1.0001*np.exp(ys)-1.0001)/(1+np.exp(ys))


class ForbesAnalysis(Analysis):
    # ylabel = "log(score)"
    title = "Forbes"

    def analyze(self):
        model = WLS(self.ys, sm.add_constant(self.xs),
                    weights=self.xs[:, -1])
        results = model.fit()
        results.coef_ = results.params[1:]
        results.intercept_ = results.params[0]
        self.model = results

    def _transform_y(self, ys):
        self.ylabel = "log(score)"
        return np.log(ys+1)

    def _inverse_y(self, ys):
        return np.exp(ys)-1

    def size_measure(self, size1, size2):
        self.xlabel = "size1+size2"
        assert size1 >= 0 and size2 >= 0
        return size1+size2


class TetraAnalysis(Analysis):
    # xlabel = "size1 + size2"
    # ylabel = "atanh(score)"
    title = "Tetrachoric correlation"

    # def _transform_y(self, ys):
    #     return np.arctanh(ys)

    def _inverse_y(self, ys):
        return np.tanh(ys)

    def size_measure(self, size1, size2):
        self.xlabel = "size1+size2"
        assert size1 >= 0 and size2 >= 0
        return size1+size2


class TrackInfo:
    def __init__(self, filename):
        self._info = pandas.read_csv(filename, sep="\t", header=4, index_col=1)
        self._info.index = [name.split("(")[0].strip() for name in self._info.index]
        self._sizes = self._info["bpcovered"]

    def get_bp_covered(self, title):
        title = title.split("(")[0].strip()
        if title not in self._sizes.index:
            print("%s not in %s" % (title, self._sizes.index))
            raise
        return self._sizes[title]


class ScoreTable:
    def __init__(self, filename):
        self._table = pandas.read_csv(
            filename, sep="\t", header=2, index_col=0)

        self.tf_types = self._parse_names()

    def _tf_from_id(self, name):
        parts = name.split()[0].split("_")
        return parts[-1]

    def _parse_names(self):
        return [self._tf_from_id(name) for name in self._table.columns]

    @classmethod
    def from_file(cls, file_name):
        with open(file_name) as f:
            first_line = f.readline()
            header = None
            while (not first_line.strip()) or first_line.startswith("#"):
                header = first_line
                first_line = f.readline()
            lines = [[float(v) for v in first_line.split()[1:]]]
            lines.extend([[float(v) for v in line.split()[1:]] for line in f])
        return cls(header, lines)


class TetraSinglePairAnalysis(Analysis):

    def _transform_y(self, ys):
        return np.arctanh(ys)

    def size_measure(self, size1, size2):
        assert size1 >= 0 and size2 >= 0
        return size1+size2

    def endcode(self):
        pass

    def analyze(self):
        pair_ids = np.array(self.pair_ids)
        pair_keys = np.unique(pair_ids)
        sizes = np.array(self.size_measures)
        scores = np.array(self.scores)
        self.coefs = []
        for pair_key in pair_keys:
            model = linear_model.LinearRegression()
            model.fit(sizes[pair_ids == pair_key][:, None],
                      scores[pair_ids == pair_key][:, None])
            print(model.coef_)
            self.coefs.append(model.coef_[0])

    def plot(self):
        ys = np.array(self.coefs)
        plt.hist(ys/self.sd)
        np.savetxt("hist.png", ys/self.sd, delimiter="\t")
        # plt.xlim(-10000, 10000)
        # plt.ylim(np.min(ys), np.max(ys))

        # plt.xlim(-102000, 10000)
        # plt.ylim(0, 1)
        plt.savefig("slopes.png")


class JaccardSinglePairAnalysis(TetraSinglePairAnalysis):
    def _size_measure(self, size1, size2):
        return np.log(size1)+np.log(size2)

    def size_measure(self, size1, size2):
        assert size1 >= 0 and size2 >= 0
        return np.log(size1*size2+1)


class ForbesSinglePairAnalysis(TetraSinglePairAnalysis):

    def _transform_y(self, ys):
        return np.log(ys+1)

    def size_measure(self, size1, size2):
        assert size1 >= 0 and size2 >= 0
        return size1+size2

    def analyze(self):
        pair_ids = np.array(self.pair_ids)
        pair_keys = np.unique(pair_ids)
        sizes = np.array(self.size_measures)
        scores = np.array(self.scores)
        self.coefs = []
        for pair_key in pair_keys:
            mask = pair_ids == pair_key
            model = WLS(scores[mask], sm.add_constant(sizes[mask]),
                        weights=sizes[mask])
            res = model.fit()
            self.coefs.append(res.params[1])


class_dict = {'tetra': [TetraAnalysis, TetraSinglePairAnalysis],
              'jaccard': [JaccardAnalysis, JaccardSinglePairAnalysis],
              'forbes': [ForbesAnalysis, ForbesSinglePairAnalysis]}


if __name__ == "__main__":
    name = 'jaccard'
    score_table = ScoreTable("scores.txt")
    track_info = TrackInfo(
        "../Galaxy52-[GSuite_(29_-_Filter_categorical_gSuite)_-_ready_for_analysis].gsuite")
    for cls in class_dict[name]:
        analysis = cls(score_table, track_info)
        analysis.generate_data()
        analysis.analyze()
        analysis.plot()
