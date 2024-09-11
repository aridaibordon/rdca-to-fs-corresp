import os

from fac_utils.compute import compute_lev

from src.utils import ELEM, ELEM_Z, default_fpath


def read_ions_elec_config(n_elec):
    with open(f"data/ions/ion{n_elec}.dat") as f:
        elec_config = [line.rstrip("\n") for line in f]
    return elec_config


def generate_atomic_data(n_elec: int, elem: str = ELEM, elem_z: int = ELEM_Z) -> None:
    path = "data/atomicdata"
    e_config = read_ions_elec_config(n_elec)

    ion = elem_z - n_elec

    fname_rdca_lev, fname_rdca_tr = default_fpath(elem, ion, "RDCA")
    fname_fs_lev, fname_fs_tr = default_fpath(elem, ion, "FS")

    compute_lev(e_config, elem, "FS", 0, path)
    compute_lev(e_config, elem, "RDCA", 0, path)

    os.system(f"mv {path}/{elem}_FS_0.lev {fname_fs_lev}")
    os.system(f"mv {path}/{elem}_RDCA_0.lev {fname_rdca_lev}")
    os.system(f"mv {path}/{elem}_FS_0.tr {path}/{fname_fs_tr}")
    os.system(f"mv {path}/{elem}_RDCA_0.tr {fname_rdca_tr}")
    os.system(f"rm {path}/*.b")
