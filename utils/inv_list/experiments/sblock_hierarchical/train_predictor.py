"""Fit integer wave selectors and report out-of-fold regret in padded bytes."""
import hashlib
import sys
from pathlib import Path

NAMES = "n blocks runs span peak gap_blocks gap_words bit_changes gap0 gap1 gap2 gap3".split()


def fold(ordinal):
    """Keep paired vectors and all their superblocks in the same fold."""
    return hashlib.blake2s(str(ordinal).encode(), digest_size=4).digest()[0] % 5


def load(directory):
    rows = []
    for dataset in ("sorted", "unsorted"):
        with (directory / (dataset + ".features")).open() as stream:
            assert next(stream).split()[2:-3] == NAMES
            for line in stream:
                ordinal, sb, *features, rival, flat, hier = map(int, line.split())
                rows.append((tuple(v.bit_length() for v in features), flat-hier,
                             fold(ordinal), dataset, features[1]))
    return rows


def fit_tree(rows, depth):
    """Choose a split by the decrease in total byte regret of its two leaves."""
    total = sum(y for x, y in rows)
    node = int(total > 0)
    if not depth or len(rows) < 64:
        return node
    best = abs(total)
    split = None
    for j in range(len(NAMES)):
        bins = {}
        for x, y in rows:
            count, weight = bins.get(x[j], (0, 0))
            bins[x[j]] = count + 1, weight + y
        count = weight = 0
        for threshold, (n, s) in sorted(bins.items()):
            count += n
            weight += s
            gain = abs(weight) + abs(total-weight)
            if min(count, len(rows)-count) >= 32 and gain > best:
                best, split = gain, (j, threshold)
    if split is None:
        return node
    j, threshold = split
    return (j, threshold,
            fit_tree([(x,y) for x,y in rows if x[j]<=threshold], depth-1),
            fit_tree([(x,y) for x,y in rows if x[j]>threshold], depth-1))


def predict_tree(node, x):
    while isinstance(node, tuple):
        j, threshold, left, right = node
        node = left if x[j] <= threshold else right
    return node


def fit_score(rows):
    """Integer coordinate search minimizing byte regret, ties prefer 108."""
    weights = [0]*(len(NAMES)+1)
    # Rows with equal coded costs cannot influence training.
    active = [(x+(1,),y) for x,y in rows if y]
    # Start from the best one-feature threshold so that a useful
    # coefficient and intercept can enter the model together.
    best_gain = 0
    for j in range(len(NAMES)):
        bins = {}
        for x, y in active:
            bins[x[j]] = bins.get(x[j], 0) + y
        total = sum(bins.values())
        lower = 0
        for threshold in sorted(bins):
            lower += bins[threshold]
            for sign, gain in ((1, total-lower), (-1, lower)):
                if gain > best_gain:
                    best_gain = gain
                    weights = [0]*(len(NAMES)+1)
                    weights[j] = sign*8
                    weights[-1] = -sign*8*threshold + (1 if sign < 0 else 0)
    scores = [sum(a*b for a,b in zip(weights,x)) for x,y in active]
    for step in (8,4,2,1):
        for iteration in range(2):
            for j in range(len(weights)):
                best = change = 0
                for move in (-step,step):
                    gain = sum(y*(int(s+move*x[j]>0)-int(s>0))
                               for (x,y),s in zip(active,scores))
                    if gain > best:
                        best, change = gain, move
                if change:
                    weights[j] += change
                    scores = [s+change*x[j] for (x,y),s in zip(active,scores)]
    return weights


def loss(y, hierarchical):
    return max(0, -y if hierarchical else y)


def main(directory):
    rows = load(directory)
    for dataset in ("sorted","unsorted"):
        print(dataset, "groups", sum(r[3]==dataset for r in rows))
    totals = {}
    for k in range(5):
        train = [(x,y) for x,y,f,d,b in rows if f!=k]
        tree, weights = fit_tree(train,3), fit_score(train)
        for x,y,f,d,b in rows:
            if f != k: continue
            policy = dict(gate16=loss(y,True) if b<=16 else 0,
                          weighted=loss(y,sum(a*v for a,v in zip(weights,x))+weights[-1]>0),
                          tree=loss(y,predict_tree(tree,x)), always108=loss(y,False))
            for name, value in policy.items():
                totals[d,name] = totals.get((d,name),0)+value
    for key, value in sorted(totals.items()):
        print(*key, value)
    train = [(x,y) for x,y,f,d,b in rows]
    print("features (integer bit widths)", *NAMES, "intercept")
    print("weights", fit_score(train))
    print("tree (feature index, threshold, <= branch, > branch), leaf 0=108 1=109")
    print(fit_tree(train,3))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
