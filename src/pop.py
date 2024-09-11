import numpy as np

from src.config import abako_to_atomic_config
from src.corresp import get_state_id, get_corresponding_states
from src.utils import ELEM_Z, get_atomicdata_fnames

from abako_utils.read import read_lp


def get_fs_corresp_id(ion, id):
    return f"{ion}.{id}"


def are_equivalent(state1, state2):
    """Check equivalence between two given states"""
    keys = ["ncomplex", "sname", "name"]
    return all([state1[key] == state2[key] for key in keys])


def compute_fs_populations(
    abako_state, corresp, temp
) -> tuple[list[float], list[dict]]:
    pop, state = abako_state["pop"], abako_state["state"]

    # Find corresponding states
    cstate = abako_to_atomic_config(state)

    id = get_state_id(cstate)
    rdca_state, fs_corresp = corresp[id]

    if len(fs_corresp) == 0:
        return [0], []

    # Compute corrected degeneracies
    subshell_deg = np.asarray([int(state["deg"]) for state in fs_corresp])

    g_rdca = rdca_state["deg"]
    g_fs = sum(subshell_deg)

    corrected_deg = (g_rdca / g_fs) * subshell_deg
    for fs_state, cdeg in zip(fs_corresp, corrected_deg):
        fs_state["deg"] = cdeg

    # Distribute populations
    subshell_deg = np.asarray([state["deg"] for state in fs_corresp])
    subshell_energies = np.asarray([state["E"] for state in fs_corresp])
    min_energy = min(subshell_energies)

    weights = subshell_deg * np.exp(-(subshell_energies - min_energy) / temp)

    fs_pop = pop * weights / sum(weights)

    return fs_pop, fs_corresp


def get_pop_corresp(path: str, temp: float):
    abako_df = read_lp(path)

    cion = -1
    fs_corresp = []
    for _, abako_state in abako_df.iterrows():
        if ELEM_Z == abako_state["ion"]:
            fs_corresp.append(
                {
                    "abako_id": abako_state["id"],
                    "ion": ELEM_Z,
                    "rdca_id": 0,
                    "fs_id": 0,
                    "E": 0.0,
                    "deg": 1,
                    "pop": abako_state["pop"],
                }
            )

            break

        if not abako_state["ion"] == cion:
            cion = abako_state["ion"]
            rdca_fname, _, fs_fname, _ = get_atomicdata_fnames(cion)

            corresp = get_corresponding_states(rdca_fname, fs_fname)

        new_pop, fs_states = compute_fs_populations(abako_state, corresp, temp)
        for pop, state in zip(new_pop, fs_states):
            fs_corresp.append(
                {
                    "abako_id": abako_state["id"],
                    "ion": cion,
                    "rdca_id": state["rdca_id"],
                    "fs_id": state["id"],
                    "E": state["E"],
                    "deg": state["deg"],
                    "pop": pop,
                }
            )

    return fs_corresp
