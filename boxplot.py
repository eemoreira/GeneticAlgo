import argparse
import matplotlib.pyplot as plt
import numpy as np


def read_numbers(filename):
    with open(filename, "r") as f:
        return [float(line.strip()) for line in f if line.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Generate a horizontal boxplot from a file of numbers."
    )

    parser.add_argument(
        "input_file",
        help="File containing one number per line",
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output image file (e.g. output.png)",
    )

    args = parser.parse_args()

    values = read_numbers(args.input_file)

    median = np.median(values)

    fig, ax = plt.subplots(figsize=(10, 2))

    ax.boxplot(
        values,
        vert=False,
        patch_artist=True,
    )

    ax.axvline(
        median,
        linestyle="--",
        label=f"Median = {median:.4f}",
    )

    ax.set_title("Box Plot")
    ax.set_xlabel("Values")
    ax.legend()

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, bbox_inches="tight")
        print(f"Saved plot to: {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
