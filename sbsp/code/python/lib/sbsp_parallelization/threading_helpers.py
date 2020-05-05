import logging
import threading
import time
from copy import copy
from typing import *

from sbsp_general import Environment
from sbsp_general.general import get_value

log = logging.getLogger(__name__)


class GenericThread (threading.Thread):
    def __init__(self, env, thread_id, data, *args):
        # type: (Environment, Any, Any, PipelineSBSPOptions, threading.Lock) -> None
        threading.Thread.__init__(self)
        self.env = env
        self.thread_id = thread_id
        self.gi = gi
        self.po = po
        self.lock = lock


    def run(self):
        sbsp_on_gi(self.gi, self.po)



def wait_for_all(active_threads):
    # type: (List[threading.Thread]) -> List[threading.Thread]

    for t in active_threads:
        t.join()
    return copy(active_threads)

def wait_for_any(active_threads, **kwargs):
    # type: (List[threading.Thread], Dict[str, Any]) -> Union[threading.Thread, None]

    sleep = get_value(kwargs, "sleep", 5)

    while True:
        if len(active_threads) == 0:
            return None
        for i, t in enumerate(active_threads):
            if not t.is_alive():
                del active_threads[i]
                return t

        time.sleep(sleep)



def run_one_per_thread(data, func_initializer=None, func_exec, **kwargs):
    # type: (Iterable[Any], Callable, Callable, Dict[str, Any]) -> None


    lock = threading.Lock()
    active_threads = list()
    thread_id = 0
    for gi in gil:
        logger.info("Scheduling: {}".format(gi.name))
        pd_work = os_join(env["pd-work"], gi.name, args.dn_run)
        curr_env = env.duplicate({"pd-work": pd_work})

        pf_output = os_join(pd_work, "output.csv")

        try:
            pf_t_db = clade_to_pf_db[gi.attributes["ancestor"]]
        except KeyError:
            raise ValueError("Unknown clade {}".format(gi.attributes["ancestor"]))

        po = PipelineSBSPOptions(
            curr_env, **vars(args), pf_t_db=pf_t_db, pf_output=pf_output, sbsp_options=sbsp_options,
            pbs_options=pbs_options,
        )

        # create working dir

        pf_list = os_join(pd_work, "query.list")
        mkdir_p(pd_work)

        # write genome to local list file
        GenomeInfoList([gi]).to_file(pf_list)

        # update custom options to local gi

        po['pf-q-list'] = pf_list

        thread = TransferThread(env, thread_id, gi, po, lock)
        thread.start()
        thread_id += 1

        active_threads.append(thread)

        # wait until number of active threads is low
        if len(active_threads) >= args.simultaneous_genomes:
            wait_for_any(active_threads)

        time.sleep(5)

    wait_for_all(active_threads)