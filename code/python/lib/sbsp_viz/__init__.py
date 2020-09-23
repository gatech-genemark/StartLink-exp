import os

runningOnPycharm = "PYCHARM_HOSTED" in os.environ
# os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import matplotlib
SMALL_SIZE = 16
MEDIUM_SIZE = 22
BIGGER_SIZE = 24

# matplotlib.use("pgf")
# print(matplotlib.get_cachedir())
matplotlib.rcParams.update({
    # "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
   'text.usetex': True,
    'pgf.rcfonts': False,
    # 'font.size': SMALL_SIZE,          # controls default text sizes
    # 'axes.titlesize': SMALL_SIZE,     # fontsize of the axes title
    # 'axes.labelsize': MEDIUM_SIZE,    # fontsize of the x and y labels
    # 'xtick.labelsize': SMALL_SIZE,    # fontsize of the tick labels
    # 'ytick.labelsize': SMALL_SIZE,    # fontsize of the tick labels
    # 'legend.fontsize': SMALL_SIZE,    # legend fontsize
    # 'figure.titlesize': BIGGER_SIZE,  # fontsize of the figure title
})

if not runningOnPycharm:
    matplotlib.use('Agg')

