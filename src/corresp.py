import pandas as pd

from fac_utils.read import read_lev


DEFAULT_DATA_PATH = "data/atomicdata"


def state_dict(state: pd.DataFrame) -> dict:
    return {
        "id": state["ID"],
        "E": state["E"],
        "ncomplex": state["ncomplex"],
        "sname": state["sname"],
        "name": state["name"],
        "deg": state["2J"] + 1,
    }


def get_state_id(state: dict) -> str:
    keys = ["ncomplex", "sname", "name"]
    return "-".join([state[key] for key in keys])


def fs_to_rdca_name_corresp(fs_level: str) -> str:
    """Return the corresponding RDCA energy level name for
    a given fine structure energy level name.
    """
    return ".".join([elem.split("(")[0] for elem in fs_level.split(".")])


def is_correspondant(fs_state: dict, rdca_state: dict) -> bool:
    """Check if a pair RDCA to FS are is_correspondant."""
    for attr in ["ncomplex", "sname"]:
        if not fs_state[attr] == rdca_state[attr]:
            return False

    if not fs_to_rdca_name_corresp(fs_state["name"]) == rdca_state["name"]:
        return False

    return True


def get_shell_mask(rdca_state: dict, fs_df: pd.DataFrame) -> list[bool]:
    """Return mask for subshell for a given RDCA state."""
    mask_ncomplex = fs_df["ncomplex"] == rdca_state["ncomplex"]
    mask_sname = fs_df["sname"] == rdca_state["sname"]

    return [bool1 and bool2 for bool1, bool2 in zip(mask_ncomplex, mask_sname)]


def get_corresponding_states(
    path_to_rdca_file: str, path_to_fs_file: str
) -> list[dict, list]:
    """Return a list that containts a dictionary with the
    information regarding the RDCA state and a list with
    all its corresponding FS states in dict format.
    """
    fs_df = read_lev(path_to_fs_file)
    rdca_df = read_lev(path_to_rdca_file)

    corresp = []
    for _, rdca_state in rdca_df.iterrows():
        shell_mask = get_shell_mask(rdca_state, fs_df)

        corresp_levels = []
        for _, fs_state in fs_df[shell_mask].iterrows():
            if is_correspondant(fs_state, rdca_state):
                fs_dict = state_dict(fs_state)
                fs_dict["rdca_id"] = rdca_state["ID"]
                corresp_levels.append(fs_dict)

        corresp.append((state_dict(rdca_state), corresp_levels))

    return {
        f"{get_state_id(rdca_state)}": (rdca_state, fs_corresp)
        for rdca_state, fs_corresp in corresp
    }
