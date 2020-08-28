# Author: Karl Gemayel
#
# Description: Save and load objects

import pickle


def save_obj(obj, name):
    # type: (object, str) -> None
    if not name.endswith(".pkl"):
        name += ".pkl"
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    # type: (str) -> object
    if not name.endswith(".pkl"):
        name += ".pkl"
    with open(name, 'rb') as f:
        loaded_data = pickle.load(f)
        return loaded_data

