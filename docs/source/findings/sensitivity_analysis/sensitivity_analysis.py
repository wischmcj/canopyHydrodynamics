from __future__ import annotations

import functools
import os
import sys
from itertools import product
from time import time

sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())
# print(sys.path)


import logging

from canopyhydro.configuration import test_input_dir
from canopyhydro.CylinderCollection import (pickle_collection,
                                            unpickle_collection)
from canopyhydro.Forester import Forester

already_run = []

log = logging.getLogger("script")


def log(_func=None, *, my_logger: logging.Logger = None):
    def decorator_log(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            args_repr = [repr(a) for a in args]
            kwargs_repr = [f"{k}={v!r}" for k, v in kwargs.items()]
            signature = ", ".join(args_repr + kwargs_repr)
            log.debug(f"function {func.__name__} called with args {signature}")
            try:
                result = func(*args, **kwargs)
                return result
            except Exception as e:
                log.exception(
                    f"Exception raised in {func.__name__}. exception: {str(e)}"
                )
                raise e

        return wrapper

    if _func is None:
        return decorator_log
    else:
        return decorator_log(_func)


@log
def initialize_collection(file="5_SmallTree", from_pickle=False, **kwargs):
    if from_pickle:
        basic_collection = load_from_pickle(**kwargs)
    else:
        log.info("initializing collection...")
        forest = Forester(test_input_dir)
        forest.get_file_names()
        forest.qsm_to_collection(file_name=file)
        basic_collection = forest.cylinder_collections[0]
        basic_collection.project_cylinders("XY")
        log.info("successfully initialized collection")
        pickle_collection(basic_collection, file)
    return basic_collection


@log
def prep_for_stats(collection, case_angle, case_name, calculate: bool = True):
    collection.initialize_digraph_from(in_flow_grade_lim=case_angle)
    collection.find_flow_components()
    if calculate:
        collection.calculate_flows()
    log.info(f"successfully prepped for stats {case_name}")


@log
def generate_statistics(collection, case_name):
    log.info(
        f"attempting to generate stats for file {collection.file_name}, case_name {case_name}"
    )
    statistics = collection.statistics(file_name_ext=case_name)
    log.info(f"attempting to generate flow file for {case_name}")
    collection.generate_flow_file(file_name_ext=case_name)


@log
def run_test_cases(
    cases_to_run, stats: bool = True, fig: bool = False, from_pickle: bool = False
):
    start = time()
    collection = None
    case_name = "inital_case_name"
    ret = []
    for case in cases_to_run:
        file_name, angle = case
        case_name = f"{angle}"
        if not collection:
            log.info(f"initializing collection for {file_name}")
            collection = initialize_collection(file_name)
            if not collection:
                dur = time() - start
                ret.append((None, f"{file_name}_{case_name}", dur))
                continue
        log.info(f"running case {file_name}_{case_name}")

        preped = prep_for_stats(collection, angle, case_name, calculate=stats)
        if stats:
            if preped:
                generate_statistics(collection, f"{case_name}")
            else:
                log.info("Error prepping, pickling and ending ")
                pickle_collection(collection, f"_prep_{case_name}")
                dur = time() - start
                ret.append((None, f"{file_name}_{case_name}", dur))
                continue
        if fig:
            draw_case(collection, angle=angle, file=file_name)
        try:
            pickle_collection(collection, f"_stats_{case_name}")
        except Exception as e:
            log.info(f"Error pickling collection for case {case_name}: {e}")
        collection = None
        dur = time() - start
        ret.append((True, f"{file_name}_{case_name}", dur))
    return ret


@log
def load_from_pickle(file, pickle_point, angle):
    pickle_file = f"{file}_pickle__{pickle_point}_{angle}"
    collection = None
    log.info(
        f"Loading collection for {file}, case {angle} from pickle file {pickle_file}"
    )
    collection = unpickle_collection(pickle_file)
    return collection


@log
def draw_case(collection=None, file: str = "", pickle_point: str = "", angle=""):
    if not collection:
        collection = load_from_pickle(file, pickle_point, angle)

    stem_flow_fig = collection.draw(
        plane="XZ",
        highlight_lambda=lambda: is_stem,
        save=True,
        file_name_ext=f"{file}_{angle}_XZ.png",
        show=False,
    )
    stemflow_and_trunk_fig = collection.draw(
        plane="XZ",
        filter_lambda=lambda: is_stem,
        highlight_lambda=lambda: branch_order == 0,
        save=True,
        file_name_ext=f"{file}_{angle}_XZ.png",
        show=False,
    )


run_cases = [
    ("Secrest32-14_000000", -0.56),
    ("Secrest10-02_000000", 0.3),
    ("Secrest26-03_000000", -0.24),
    ("Secrest03-12_000000", -0.04),
    ("Secrest24-07_000000", -0.88),
    ("Secrest18-13_000000", 0.08),
]

angles = [0.32, 0.24, 0.16, 0.08, 0, -0.08, -0.16, -0.24, -0.32]


@log
def get_cases(file_names, already_run, angles_to_tests):
    print("getting_cases")
    already_run = [(x, float(y)) for x, y in already_run]
    cases = product(file_names, angles_to_tests)
    return [case for case in cases if case not in already_run]


def sensitivity_analysis():
    files_to_test = ["5_SmallTree"]
    cases_to_run = get_cases(files_to_test, run_cases, angles)
    log.info(f"Will run {len(cases_to_run)} cases : {cases_to_run}")
    print(f"running cases {cases_to_run}")
    run_test_cases(cases_to_run, fig=False)


if __name__ == "__main__":
    log = logging.getLogger("model")
    # load_from_pickle('5_SmallTree_1', 'stats', 0.36666)
    log.info("starting")
    sensitivity_analysis()
