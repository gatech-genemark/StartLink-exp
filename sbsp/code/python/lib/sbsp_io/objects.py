# Author: Karl Gemayel
#
# Description: Save and load objects

import pickle


def save_obj(obj, name):
    # type: (object, str) -> None
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    # type: (str) -> object
    with open(name + '.pkl', 'rb') as f:
        loaded_data = pickle.load(f)
        return loaded_data

