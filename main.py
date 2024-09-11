import os
import re

from src.transdata import create_transdata_file


RDCA_PATH = "/Banco/aridai/rdca_to_fs/out/rdca"
FS_PATH = "/Banco/aridai/rdca_to_fs/out/fs"


def is_transdata_file(fname: str) -> bool:
    return bool(re.search("transdata", fname))


def is_radiative_file(fname: str) -> bool:
    return bool(re.search("em_op", fname))


def main() -> None:
    print("Computing correspondence mechanism...")

    files = os.listdir(RDCA_PATH)
    files.sort()

    for fname in files:
        if not is_transdata_file(fname):
            continue
        
        path_fs_trans = os.path.join(FS_PATH, fname)
        path_rdca_trans = os.path.join(RDCA_PATH, fname)
        path_rdca_lp = path_rdca_trans.replace("transdata_Ar", "Ar_lp_crm")
        
        create_transdata_file(path_fs_trans, path_rdca_lp, path_rdca_trans)

        print(f"Generated FS for {fname}")


if __name__ == "__main__":
    main()
