import os
import pandas as pd
import glob
from biopandas.pdb import PandasPdb
from tqdm import tqdm


def main() -> None:
    os.chdir("pdb_data")
    exist_pdb_files = glob.glob("*.pdb")
    pdb_pair_df = pd.read_csv("pdb_with_hinges.csv")
    pdb_pair_df["p_pdb_id"] = pdb_pair_df["p_pdb"].apply(lambda x: x.split("_")[0]).str.lower()
    pdb_pair_df["q_pdb_id"] = pdb_pair_df["q_pdb"].apply(lambda x: x.split("_")[0]).str.lower()
    pdb_pair_df["p_chain_id"] = pdb_pair_df["p_pdb"].apply(lambda x: x.split("_")[1])
    pdb_pair_df["q_chain_id"] = pdb_pair_df["q_pdb"].apply(lambda x: x.split("_")[1])
    ok_list = []
    ng_list = []
    for pdb_1, pdb_2, chain_1, chain_2 in pdb_pair_df[["p_pdb_id", "q_pdb_id", "p_chain_id", "q_chain_id"]].to_numpy():
        diff_len = load_and_output_pdb(
            pdb_id_1=pdb_1, chain_id_1=chain_1, pdb_id_2=pdb_2, chain_id_2=chain_2, exist_pdb_files=exist_pdb_files
        )
        print(diff_len)
        if diff_len > 0:
            ng_list.append((pdb_1, pdb_2, chain_1, chain_2, diff_len))
        else:
            ok_list.append((pdb_1, pdb_2, chain_1, chain_2))
    return ok_list, ng_list

def load_df(pdb_id: str, chain_id: str, exist_pdb_files):
    if not f"pdb{pdb_id}_{chain_id}.pdb" in exist_pdb_files:
        pdb_df = PandasPdb().fetch_pdb(pdb_id)
        pdb_df.to_pdb(path=f"pdb{pdb_id}_{chain_id}.pdb", records=None, gz=False, append_newline=True)
        df = pdb_df.df["ATOM"][
            [
                "atom_name",
                "residue_number",
                "residue_name",
                "chain_id",
                "x_coord",
                "y_coord",
                "z_coord",
            ]
        ]
    else:
        df = (
            PandasPdb()
            .read_pdb(path=f"pdb{pdb_id}_{chain_id}.pdb")
            .df["ATOM"][
                [
                    "atom_name",
                    "residue_number",
                    "residue_name",
                    "chain_id",
                    "x_coord",
                    "y_coord",
                    "z_coord",
                ]
            ]
        )
    if chain_id != "":
        df = df[df["chain_id"] == chain_id]
    else:
        df = df[df["chain_id"] == df["chain_id"].unique()[0]]
    df = df[df["atom_name"] == "CA"]
    df = df.drop_duplicates(subset=["residue_number"])
    return df


def load_and_output_pdb(pdb_id_1: str, chain_id_1: str, pdb_id_2: str, chain_id_2: str, exist_pdb_files):
    df_1 = load_df(pdb_id=pdb_id_1, chain_id=chain_id_1, exist_pdb_files=exist_pdb_files)
    df_2 = load_df(pdb_id=pdb_id_2, chain_id=chain_id_2, exist_pdb_files=exist_pdb_files)
    df_1_merge = df_1.merge(df_2["residue_number"], on=["residue_number"], how="inner")
    df_2_merge = df_2.merge(df_1["residue_number"], on=["residue_number"], how="inner")
    df_1_merge[["x_coord", "y_coord", "z_coord"]].T.to_csv(
        f"../coord_csv/coord_{pdb_id_1}_{chain_id_1}_{pdb_id_2}_{chain_id_2}.csv",
        index=False,
        header=None,
    )
    df_2_merge[["x_coord", "y_coord", "z_coord"]].T.to_csv(
        f"../coord_csv/coord_{pdb_id_2}_{chain_id_2}_{pdb_id_1}_{chain_id_1}.csv",
        index=False,
        header=None,
    )
    is_aligned = len(df_1_merge[df_1_merge["residue_name"] == df_2_merge["residue_name"]]) == len(df_1_merge)
    if not is_aligned:
        print(df_1_merge[df_1_merge["residue_name"] != df_2_merge["residue_name"]]["residue_name"])
        print(df_2_merge[df_1_merge["residue_name"] != df_2_merge["residue_name"]]["residue_name"])
    if is_aligned:
        diff_len = 0
    else:
        diff_len = len(df_1_merge[df_1_merge["residue_name"] != df_2_merge["residue_name"]])
    return diff_len

if __name__ == "__main__":
    ok, ng = main()
    pd.DataFrame(ok, columns=["p_pdb_id", "q_pdb_id", "p_chain_id", "q_chain_id"]).to_csv("ok_pdb.csv", index=False)
