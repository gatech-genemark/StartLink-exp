# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 2/24/20

import logging
import argparse
import os
import threading
import random
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from sbsp_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from sbsp_general.general import run_shell_cmd
from sbsp_io.general import write_string_to_file

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-data', required=True, help="File to be transfered")
parser.add_argument('--pd-compute-destination', required=True, help="Location where file will be transfered on compute nodes")

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


def create_pbs_file_for_transfer(env, pf_pbs, sender, receiver, pf_data_from, pf_data_to):
    # type: (Environment, str, str, str, str, str) -> None

    pd_data_to = os.path.dirname(pf_data_to)

    cmd = ""
    cmd += "#PBS -N {}\n".format(receiver)
    cmd += "#PBS -o {}\n".format(os.path.join(env["pd-work"], "{}.oe".format(receiver)))
    cmd += "#PBS -j oe\n"
    cmd += "#PBS -l nodes={}\n".format(receiver)
    cmd += "#PBS -l walltime=07:00:00:00\n"
    cmd += "#PBS -V\n"
    cmd += "#PBS -W umask=002\n"

    cmd += "mkdir -p {}\n".format(pd_data_to)
    cmd += "PBS_O_WORKDIR={}\n".format(pd_data_to)
    cmd += "cd $PBS_O_WORKDIR\n"

    if sender is None:
        cmd += "rsync {} {}\n".format(sender, pf_data_from, pf_data_to)
    else:
        cmd += "rsync {}:{} {}\n".format(sender, pf_data_from, pf_data_to)

    write_string_to_file(cmd, pf_pbs)




def transfer_from_compute_to_compute(env, sender, receiver, pf_data_from, pf_data_to, fn_tmp_prefix=None):
    # type: (Environment, str, str, str, str, str) -> None

    if fn_tmp_prefix is None:
        fn_tmp_prefix = ""

    pf_pbs = os.path.join(env["pd-work"], "{}_transfer.pbs".format(fn_tmp_prefix))
    create_pbs_file_for_transfer(env, pf_pbs, sender, receiver, pf_data_from, pf_data_to)
    return
    job_id = run_shell_cmd("qsub -V {}".format(pf_pbs))

    cmd = r'while [ $(qstat -a | grep " R\|Q\|H " | grep ' + job_id + r'  | wc -l) != 0 ]; do sleep 60 ; done'

    # wait for job to finish
    run_shell_cmd(cmd, do_not_log=True)

def transfer_from_head_to_compute(env, receiver, pf_data_from, pf_data_to, fn_tmp_prefix=None):
    # type: (Environment, str, str, str, str) -> None

    if fn_tmp_prefix is None:
        fn_tmp_prefix = ""

    pf_pbs = os.path.join(env["pd-work"], "{}_transfer.pbs".format(fn_tmp_prefix))
    create_pbs_file_for_transfer(env, pf_pbs, None, receiver, pf_data_from, pf_data_to)
    return
    job_id = run_shell_cmd("qsub -V {}".format(pf_pbs))

    cmd = r'while [ $(qstat -a | grep " R\|Q\|H " | grep ' + job_id + r'  | wc -l) != 0 ]; do sleep 60 ; done'

    # wait for job to finish
    run_shell_cmd(cmd, do_not_log=True)



class TransferThread (threading.Thread):
    def __init__(self, env, thread_id, pf_data, sender, receiver, lock):
        # type: (Environment, Any, str, str, str, threading.Lock) -> None
        threading.Thread.__init__(self)
        self.env = env
        self.thread_id = thread_id
        self.pf_data = pf_data
        self.sender = sender
        self.receiver = receiver
        self.lock = lock


    def run(self):
        transfer_from_compute_to_compute(
            self.env, self.sender, self.receiver, self.pf_data, self.pf_data, fn_tmp_prefix=self.thread_id
        )

        run_shell_cmd("sleep {}".format(random.randint(1, 2)))


def get_up_nodes():
    # type: () -> List[str]
    return ["node{}".format(x) for x in range(20)]
    output = run_shell_cmd("pbsnodes -l up | awk '{print $1}'", True)
    return output.strip().split("\n")


def smart_transfer(env, pf_data, pd_destination):
    # type: (Environment, str, str) -> None

    list_available_nodes = get_up_nodes()

    list_can_send = list()
    list_can_receive = list_available_nodes.copy()

    # first transfer is from head node to a compute node
    receiver = list_can_receive.pop()
    pf_compute_data = os.path.join(pd_destination, os.path.basename(pf_data))


    transfer_from_head_to_compute(env, receiver, pf_data, pf_compute_data, "head")
    list_can_send.append(receiver)

    lock = threading.Lock()

    active_threads = set()
    thread_id = 0

    order_of_completion = list()
    # run multithreaded, copying
    while len(list_can_receive) > 0:

        while len(list_can_send) > 0 and len(list_can_receive) > 0:
            receiver = list_can_receive.pop()
            sender = list_can_send.pop()

            if sender is None:
                raise RuntimeError("No sender available, something went wrong.")

            thread = TransferThread(env, thread_id, pf_compute_data, sender, receiver, lock)
            thread.start()


            active_threads.add(thread)

        thread_id += 1

        # wait until at least one thread is done
        completed_threads = set()
        while True:
            for t in active_threads:
                if not t.is_alive():
                    completed_threads.add(t)

                    # add to lists
                    list_can_send.append(t.receiver)
                    list_can_send.append(t.sender)

                    order_of_completion.append((t.sender, t.receiver, t.thread_id))

            if len(completed_threads) > 0:
                for t in completed_threads:
                    active_threads.remove(t)
                break

    for t in active_threads:
        t.join()
        order_of_completion.append((t.sender, t.receiver, t.thread_id))

    for o in sorted(order_of_completion, key=lambda x: x[2]):
        print("{}: {} -> {}".format(o[2], o[0], o[1]))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    smart_transfer(env, os.path.abspath(args.pf_data), os.path.abspath(args.pd_compute_destination))
    pass


if __name__ == "__main__":
    main(my_env, parsed_args)
