import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse


def read_ga_file(filename):
    metadata = {}
    data_start = 0

    with open(filename, "r") as f:
        lines = f.readlines()

    # extract metadata
    for i, line in enumerate(lines):
        line = line.strip()

        if line.startswith("Generations:"):
            metadata["GEN"] = line.split(":")[1].strip()
        elif line.startswith("Population Size:"):
            metadata["POP"] = line.split(":")[1].strip()
        elif line.startswith("Mutation Rate:"):
            metadata["PM"] = line.split(":")[1].strip()
        elif line.startswith("Crossover Rate:"):
            metadata["PC"] = line.split(":")[1].strip()
        elif line.startswith("Elitism Selection:"):
            metadata["ELIT"] = line.split(":")[1].strip()
        elif line.startswith("Generation"):
            data_start = i
            break

    df = pd.read_csv(
        filename,
        skiprows=data_start
    )

    return df, metadata


def plot(files, output):
    plt.figure()

    for file in files:
        df, meta = read_ga_file(file)

        x = df["Generation"].values
        best = df["BestFitness"].values
        avg = df["AverageFitness"].values

        label = (
            f"POP={meta['POP']} "
            f"PC={meta['PC']} "
            f"PM={meta['PM']} "
            f"ELIT={meta['ELIT']}"
        )

        line, = plt.plot(
            x,
            best,
            label=f"{label} (Best)"
        )
        color = line.get_color()

        plt.scatter(x, best, s=10, color=color)

        plt.plot(
            x,
            avg,
            color="red",
            label=f"{label} (Average)"
        )

        plt.scatter(
            x,
            avg,
            s=10,
            color="red"
        )

    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.title("GA Convergence Comparison")

    plt.legend()
    plt.grid()

    plt.savefig(output)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot GA convergence")

    parser.add_argument(
        "files",
        nargs="+",
        help="Input files (e.g., output.txt)"
    )

    parser.add_argument(
        "-o", "--output",
        default="convergence.png",
        help="Output image filename"
    )

    args = parser.parse_args()

    plot(args.files, args.output)
