#!/usr/bin/env python3
"""
aaKomp score calculation script
Written by Johnathan Wong
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def parse_file(fname):
    d = {}
    with open(fname) as f:
        next(f)  # skip header
        for line in f:
            cols = line.strip().split()
            if len(cols) < 9:
                continue
            try:
                val = float(cols[5])
            except ValueError:
                continue
            key = cols[8]
            d[key] = max(d.get(key, 0), val)
    return list(d.values())

def make_cdf(vals):
    x = np.sort(vals)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

def save_plot(x, y, out_png):
    plt.figure()
    plt.plot(x, y, drawstyle='steps-post')
    plt.xlabel("Value")
    plt.ylabel("CDF")
    plt.title("aaKomp CDF")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python aakomp_score.py <input_file>")
        sys.exit(1)

    in_file = sys.argv[1]
    base = os.path.basename(in_file)

    vals = parse_file(in_file)
    if not vals:
        print(f"Nothing to score in {base}")
        sys.exit(1)

    x, y = make_cdf(vals)
    auc = np.trapezoid(y, x)
    score = (1 - auc) * 100

    print(f"aaKomp score for {base}: {score:.2f}")

    with open(f"{in_file}_score.txt", "w") as out:
        out.write(f"{score:.2f}\n")

    save_plot(x, y, f"{in_file}_cdf.png")
