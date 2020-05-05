# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import threading
import random
import time
from typing import *
import numpy as np
from multiprocessing import Pool

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger

class TransferThread (threading.Thread):
    def __init__(self, env, thread_id, results):
        # type: (Environment, Any, List[Any]) -> None
        threading.Thread.__init__(self)
        self.env = env
        self.thread_id = thread_id
        self.results = results

    def run(self):
        logger.info(f"Setting seed for thread {self.thread_id}")
        random.seed(1)

        l = list()
        for r in range(10):
            l.append(random.random())
            time.sleep(0.001)

        self.results[self.thread_id] = l



def wait_for_all(active_threads):
    # type: (List[T]) -> List[T]

    done_threads = list()
    while True:
        if len(active_threads) == 0:
            break
        for i, t in enumerate(active_threads):
            if not t.is_alive():
                del active_threads[i]
                done_threads.append(t)

        time.sleep(5)

    return done_threads

def lists_are_equal(lists):
    # type:(List[List[Any]]) -> bool

    if len(lists) == 0:
        return True

    list_len = len(lists[0])

    # check all same length
    for i in range(1, len(lists)):
        if len(lists[i]) != list_len:
            return False

    for j in range(list_len):
        element = lists[0][j]

        for i in range(1, len(lists)):
            if lists[i][j] != element:
                return False

    return True



def test1(env):
    # type: (Environment) -> bool

    # Sketch:
    # Test whether random seeds are affected by each other in multiple threads
    num_threads = 20
    active_threads = list()     # type: List[TransferThread]
    values = [None] * num_threads

    for i in range(num_threads):
        active_threads.append(
            TransferThread(env, i, values)
        )
        active_threads[-1].start()

    done_threads = wait_for_all(active_threads)     # type: List[TransferThread]

    for l in values:
        print (l)
    return lists_are_equal(values)


def Foo_np(seed=None):
    x = np.random.RandomState(seed)
    # if not seed:
    #     np.random.seed(seed)
    # return np.random.uniform(0, 1, 10)
    return x.uniform(0, 1, 10)

def test2(env):
    # type: (Environment) -> bool
    pool = Pool(processes=8)
    x = np.array(pool.map(Foo_np, [1]*30))
    print(x)

def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    # if test1(env):
    #     logger.info("Test 1 successful")
    # else:
    #     logger.info("Test 1 unsuccessful")

    test2(env)


if __name__ == "__main__":
    main(my_env, parsed_args)
