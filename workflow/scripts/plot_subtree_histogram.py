# -*- coding: utf-8 -*-
import argparse
import ete3
from src.treeinform_collapse import load_trees
import matplotlib.pyplot as plt


def is_ultrametric(tree):
    """
    https://www.biostars.org/p/297625/
    """
    last_distance = None
    distance = 0
    for post, node in tree.iter_prepostorder():
        if post:
            distance -= node.dist
        else:
            if not node.is_leaf():
                distance += node.dist
            else:
                if last_distance is None:
                    print(last_distance)
                    last_distance = distance + node.dist
                elif last_distance != distance + node.dist:
                    return False, last_distance
    return True, last_distance


def branch_length_histogram(newicks):
    """
    # Agalma - an automated phylogenomics workflow
    # Copyright (c) 2012-2017 Brown University. All rights reserved.
    # Modified by Natasha Picciani on Oct 18

    Distribution of subtree branch lengths. Returns a dictionary
    containing subtree lengths and corresponding counts.

    Args:
    -- newicks: list containing the absolute paths to tree files
    """
    hist = {}
    for newick in newicks:
        tree = ete3.Tree(newick)
        outgroup = tree.get_midpoint_outgroup()
        if not outgroup is None:
            tree.set_outgroup(outgroup)
        tree.convert_to_ultrametric(
            tree_length=1
        )  # converts tree to ultrametric with length 1
        if is_ultrametric(tree) == False:
            print("tree is not ultrametric")
        for node in tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node.add_feature("branchlength", 0)
                node.add_feature("under", True)
            if not node.is_leaf():
                children = node.get_children()
                branchlength = (
                    children[0].get_distance(children[1])
                    + children[0].branchlength
                    + children[1].branchlength
                )
                node.add_feature("branchlength", branchlength)
                hist[branchlength] = hist.get(branchlength, 0) + 1
    return hist


def plot_histogram(histogram, binsize, threshold_value, outdir):
    """ "
    Plot distribution of subtree lengths.
    Returns png image with a histogram of the distribution.

    Args:
    -- histogram: dictionary containing subtree lengths and corresponding counts.
    -- outdir: output directory for image file named 'branch.length.hist.png'.
    """

    figname = f"{outdir}/branch.length.hist.png"
    plt.hist(histogram, bins=binsize, edgecolor="k")
    plt.xlabel("Branch Length")
    plt.ylabel("Frequency")
    plt.title("Distribution of Subtree Branch Lengths")
    plt.axis([0, 20, 0, 1000])
    plt.axvline(threshold_value, color="red")
    plt.savefig(figname, facecolor="white")


def main(args):

    gene_trees_folder = args.gt
    outdir_name = args.o
    threshold_value = float(args.t)
    binsize = int(args.b)

    trees = load_trees(gene_trees_folder)
    histogram = branch_length_histogram(trees)
    plot_histogram(histogram, binsize, threshold_value, outdir_name)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Plot distribution of subtree branch lengths from several gene trees"
    )

    parser.add_argument(
        "-gt",
        "-genetrees",
        type=str,
        required=True,
        help="path to folder with gene trees produced with Orthofinder",
    )
    parser.add_argument(
        "-t",
        "-threshold",
        type=float,
        required=True,
        default="0.05",
        help="threshold of subtree length for collapsing gene variants",
    )
    parser.add_argument(
        "-b",
        "-binsize",
        type=int,
        required=True,
        default="10000",
        help="bin size for histogram",
    )

    parser.add_argument(
        "-o", "-outdir", type=str, required=False, default=".", help="output directory"
    )

    args = parser.parse_args()
    main(args)
