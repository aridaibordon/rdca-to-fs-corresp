import pandas as pd

from itertools import product

from src.pop import get_pop_corresp

from fac_utils.read import read_tr
from abako_utils.read import read_transdata, read_transdata_header, read_transdata_cond


def get_fs_transitions(low, upp, corresp_df) -> list[tuple[int, dict]]:
    low_df = corresp_df[corresp_df["abako_id"] == low]
    upp_df = corresp_df[corresp_df["abako_id"] == upp]

    low_fs = [(index, row.to_dict()) for index, row in low_df.iterrows()]
    upp_fs = [(index, row.to_dict()) for index, row in upp_df.iterrows()]

    return list(product(low_fs, upp_fs))


def compute_einstein_coef(
    delta_E_rdca,
    delta_E_fs,
    ldeg,
    udeg,
    ldeg_rdca,
    udeg_rdca,
    ein_rdca,
) -> float:
    return (
        (delta_E_fs / delta_E_rdca) ** 3
        * (udeg / udeg_rdca)
        * (ldeg_rdca / ldeg)
        * ein_rdca
    )


def get_einstein_coef(lind, uind, tr_df):
    tr_df = tr_df[(tr_df["upp"] == uind) & (tr_df["low"] == lind)]

    if tr_df.empty:
        return 0
    return tr_df["ein"].iloc[0]


def get_bb_transrow(lind, uind, lstate, ustate, rdca_transrow, tr_df) -> list:
    delta_E_rdca = rdca_transrow["deltaE"]
    delta_E_fs = ustate["E"] - lstate["E"]

    """
    cein = compute_einstein_coef(
        delta_E_rdca,
        delta_E_fs,
        lstate["deg"],
        ustate["deg"],
        rdca_transrow["ldeg"],
        rdca_transrow["udeg"],
        rdca_transrow["ein"],
    )
    """

    cein = get_einstein_coef(lstate["fs_id"], ustate["fs_id"], tr_df)

    return [
        rdca_transrow["ion"],
        lind,
        uind,
        lstate["deg"],
        ustate["deg"],
        lstate["pop"],
        ustate["pop"],
        delta_E_fs,
        cein,
        rdca_transrow["dop"],
        rdca_transrow["voi"],
        rdca_transrow["sta"],
    ]


def get_bf_transrow(lind, uind, lstate, ustate, rdca_transrow):
    return [
        rdca_transrow["lion"],
        rdca_transrow["uion"],
        lind,
        uind,
        lstate["deg"],
        ustate["deg"],
        lstate["pop"],
        ustate["pop"],
        rdca_transrow["unk1"],
        rdca_transrow["unk2"],
        rdca_transrow["unk3"],
        rdca_transrow["unk4"],
        rdca_transrow["unk5"],
    ]


def get_fs_transdata(path_to_lp, path_to_transdata):
    temp, _ = read_transdata_cond(path_to_transdata)
    rdca_transbb, rdca_transbf = read_transdata(path_to_transdata)
    
    corresp = get_pop_corresp(path_to_lp, temp)
    corresp_df = pd.DataFrame(corresp)

    cion = -1
    fs_transbb = []
    for _, rdca_transrow in rdca_transbb.iterrows():
        ion, low, upp = rdca_transrow["ion"], rdca_transrow["low"], rdca_transrow["upp"]

        if ion != cion:
            tr_df = read_tr(f"data/atomicdata/Ar_{18-ion}_FS_0.tr")

        fs_transitions = get_fs_transitions(low, upp, corresp_df)
        for transition in fs_transitions:
            (lind, lstate), (uind, ustate) = transition

            if ustate["E"] - lstate["E"] < 0:
                ustate, lstate = lstate, ustate

            fs_transrow = get_bb_transrow(
                lind, uind, lstate, ustate, rdca_transrow, tr_df
            )
            fs_transbb.append(fs_transrow)

    fs_transbf = []
    for _, rdca_transrow in rdca_transbf.iterrows():
        low, upp = rdca_transrow["low"], rdca_transrow["upp"]

        fs_transitions = get_fs_transitions(low, upp, corresp_df)
        for transition in fs_transitions:
            (lind, lstate), (uind, ustate) = transition

            if ustate["E"] - lstate["E"] < 0:
                ustate, lstate = lstate, ustate

            fs_transrow = get_bf_transrow(lind, uind, lstate, ustate, rdca_transrow)
            fs_transbf.append(fs_transrow)

    fkeys_bb = {
        "ion": int,
        "low": int,
        "upp": int,
        "ldeg": float,
        "udeg": float,
        "lpop": float,
        "upop": float,
        "deltaE": float,
        "ein": float,
        "dop": float,
        "voi": float,
        "sta": str,
    }
    fkeys_bf = {
        "lion": int,
        "uion": int,
        "low": int,
        "upp": int,
        "ldeg": float,
        "udeg": float,
        "lpop": float,
        "upop": float,
        "unk1": float,
        "unk2": int,
        "unk3": int,
        "unk4": int,
        "unk5": str,
    }

    return (
        pd.DataFrame(fs_transbb, columns=fkeys_bb.keys()).astype(fkeys_bb),
        pd.DataFrame(rdca_transbf, columns=fkeys_bf.keys()).astype(fkeys_bf),
    )


def create_transdata_file(
    path_fs_transdata: str, path_rdca_lp: str, path_rdca_transdata: str
) -> None:
    
    get_fs_transdata(path_rdca_lp, path_rdca_transdata)
    fs_bbtrans, fs_bftrans = get_fs_transdata(path_rdca_lp, path_rdca_transdata)

    with open(path_fs_transdata, "w") as ftrans:
        ftrans.write(f"{read_transdata_header(path_rdca_transdata)}")

        # bb transitions
        for row in fs_bbtrans.values.tolist():
            (
                ion,
                low,
                upp,
                ldeg,
                udeg,
                lpop,
                upop,
                deltaE,
                ein,
                dop,
                voi,
                sta,
            ) = row
            ftrans.write(
                f"{ion:6d} {low:6d} {upp:6d} {ldeg:6.1f} {udeg:6.1f}"
                + f" {lpop:14.6e} {upop:14.6e} {deltaE:14.6e}"
                + f" {ein:14.6e} {dop:14.6e} {voi:14.6e} {sta}\n"
            )

        ftrans.write("    -1\n\n")

        # bf transitions
        for row in fs_bftrans.values.tolist():
            (
                lion,
                uion,
                lind,
                uind,
                ldeg,
                udeg,
                lpop,
                upop,
                unk1,
                unk2,
                unk3,
                unk4,
                sta,
            ) = row

            ftrans.write(
                f"{lion:6d} {uion:6d} {lind:6d} {uind:6d}"
                + f" {ldeg:6.1f} {udeg:6.1f} {lpop:14.6e} {upop:14.6e}"
                + f" {unk1:14.6e} {unk2:3d} {unk3:3d} {unk4:8d} {sta}\n"
            )

        ftrans.write("    -1")
