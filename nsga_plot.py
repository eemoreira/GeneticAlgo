#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import os

print("CWD:", os.getcwd())

def make_plot(df, x, y, c, outdir):
    plt.figure()

    sc = plt.scatter(
        df[x],
        df[y],
        c=df[c]
    )

    plt.xlabel(x)
    plt.ylabel(y)
    plt.colorbar(sc, label=c)
    plt.tight_layout()

    outfile = outdir / f"{x}-{y}.png"
    plt.savefig(outfile)
    plt.close()

    print(f"Saved {outfile}")


def main():
    if len(sys.argv) != 3:
        print(
            f"Usage: {sys.argv[0]} <input.csv> <output_dir>",
            file=sys.stderr,
        )
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    outdir = Path(sys.argv[2])

    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)

    make_plot(df, "f1", "f2", "f3", outdir)
    make_plot(df, "f1", "f3", "f2", outdir)
    make_plot(df, "f2", "f3", "f1", outdir)


if __name__ == "__main__":
    main()
