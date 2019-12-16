import os

runningOnPycharm = "PYCHARM_HOSTED" in os.environ
if not runningOnPycharm:
    import matplotlib
    matplotlib.use('Agg')