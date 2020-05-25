import logging
from typing import *
from sbsp_general.MGMMotifModel import MGMMotifModel
import seaborn
import matplotlib.pyplot as plt

log = logging.getLogger(__name__)

class MGMMotifModelVisualizer:

    @staticmethod
    def _viz_letter(mgm_mm, letter, ax):
        # type: (MGMMotifModel, str, plt.Axes) -> None

        x = range(mgm_mm.motif_width())
        y = [mgm_mm._motif[letter][p] for p in x]

        seaborn.lineplot(x, y, ax=ax, color="blue")

    @staticmethod
    def visualize(mgm_mm, title=""):
        # type: (MGMMotifModel, str) -> None

        fig = plt.figure(figsize=(10, 12))
        shape = (4, 2)

        ax1 = plt.subplot2grid(shape, (0, 0))
        ax2 = plt.subplot2grid(shape, (0, 1))
        ax3 = plt.subplot2grid(shape, (1, 0))
        ax4 = plt.subplot2grid(shape, (1, 1))
        ax_logo = plt.subplot2grid(shape, (3, 0))
        ax_counts = plt.subplot2grid(shape, (2, 0))
        ax_pos_dist = plt.subplot2grid(shape, (2, 1))
        ax_text = plt.subplot2grid(shape, (3, 1))

        axes = [ax1, ax2, ax3, ax4]     # letters


        ylim = [-0.1, 1.1]
        letters = sorted(mgm_mm._motif.keys())

        # for each letter
        for l, ax in zip(letters, axes):
            MGMMotifModelVisualizer._viz_letter(mgm_mm, l, ax)
            ax.set_title(f"{l}")
            ax.set_ylim(ylim)

        plt.suptitle("Gc range: {}".format(title))

        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.show()