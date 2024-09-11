from os.path import join


# Default element attributes
ELEM = "Ar"
ELEM_Z = 18

# Default path
DEFAULT_ATOMICDATA_PATH = "data/atomicdata"

# Energy intervals for transdata (eV))
ENERGY_LIMITS = (3000, 5000)


def default_fpath(elem, ion, energy_scheme, ci=0) -> tuple[str]:
    return (
        join(DEFAULT_ATOMICDATA_PATH, f"{elem}_{ion}_{energy_scheme}_{ci}.lev"),
        join(DEFAULT_ATOMICDATA_PATH, f"{elem}_{ion}_{energy_scheme}_{ci}.tr"),
    )


def get_atomicdata_fnames(ion: int, elem: str = ELEM):
    """Return a tuple with the default RDCA and FS atomic datafiles
    for a given element and ion."""
    RDCA_lev, RDCA_tr = default_fpath(elem, ion, "RDCA")
    FS_lev, FS_tr = default_fpath(elem, ion, "FS")
    return RDCA_lev, RDCA_tr, FS_lev, FS_tr
