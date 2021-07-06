"""
Author: Michael Markus Ackermann
================================
"""

import pickle


def __endswith_pkl_check(file_name: str) -> None:
    """Checks if file_name endswith pickle extension.

    Args:
        file_name (str): file_name.

    Raises:
        ValueError: If the file_name doesn't ends with '.pkl'.

    Returns:
        None
    """
    if not file_name.endswith(".pkl"):
        raise ValueError(
            "file_name= {} doesn't ends with '.pkl'.".format(file_name))
    return None


def save_object_with_pickle(file_name: str, obj: object) -> None:
    """Save a object easily be using pickle module.

    Args:
        file_name (str): file_name.pkl
        obj (object): Any PyAbsorp object.

    Returns:
        None
    """
    __endswith_pkl_check(file_name)
    with open(file_name, "wb") as dump_file:
        pickle.dump(obj, dump_file, pickle.HIGHEST_PROTOCOL)
        dump_file.close()
    return None


def load_object_with_pickle(file_name: str) -> object:
    """Loads a object easily be using pickle module.

    Args:
        file_name (str): file_name.pkl

    Returns:
        object (object): Any PyAbsorp object.
    """
    __endswith_pkl_check(file_name)
    try:
        with open(file_name, "wb") as load_file:
            obj = pickle.load(load_file)
            load_file.close()
    except (pickle.UnpicklingError, ImportError, ValueError):
        print("Something went wrong in the loading process.")
        raise
    return obj
