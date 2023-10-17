import pandas as pd
import matplotlib.pyplot as plt
import os

expert_annotated_list = [
    "35",
    "25 : 55",
    "98 : 110",
    "13 : 81",
    "89 : 264",
    "135",
    "150 : 318",
    "86 : 177",
    "91 : 251",
    "91 : 193",
    "104 : 236",
    "36 : 63 : 102",
]

shibuya_method_list = [
    "35 : 52",
    "45 : 57",
    "98 : 110",
    "12 : 76",
    "84 : 253",
    "84 : 134",
    "42 : 141",
    "85 : 179",
    "92 : 251",
    "91 : 192",
    "103 : 235",
    "35 : 71",
]
flexprot_list = [
    "39",
    "",
    "",
    "",
    "89",
    "",
    "",
    "88 : 181",
    "85 : 245",
    "84 : 177",
    "101",
    "35 : 67",
]
fatcat_list = ["", "", "", "", "91", "", "", "89 : 183", "92 : 252", "81", "99", ""]
nahal_list = [
    "61",
    "48 : 59",
    "98 : 111",
    "14 : 83",
    "87 : 252",
    "72 : 136",
    "47 : 148",
    "79 : 179",
    "90 : 253",
    "90 : 162",
    "103 : 237",
    "36 : 69",
]


def main():
    os.makedirs("figures", exist_ok=True)
    small_monge_2 = pd.read_csv("rmsdhk_monge_check_result.csv")
    more_monge_2 = pd.read_csv(
        "rmsdh_result/rmsdhk_monge_check_more_data_2.csv",
        names=[
            "p_pdb_id",
            "q_pdb_id",
            "Residue length",
            "RMSD",
            "",
            "monge_rate_1",
            "monge_rate_2",
            "monotonicity_rate_1",
            "monotonicity_rate_2",
        ],
        header=None,
        skiprows=1,
    )
    small_fast_rmsdhk_2 = pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_2.csv")
    small_rmsdhk_2 = pd.read_csv("rmsdh_result/rmsdh_hingek_cnt_2.csv")
    small_fast_rmsdhk_3 = pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_3.csv")
    small_rmsdhk_3 = pd.read_csv("rmsdh_result/rmsdh_hingek_cnt_3.csv")
    small_fast_rmsdhk_4 = pd.read_csv("rmsdh_result/fast_rmsdh_hingek_cnt_4.csv")
    small_rmsdhk_4 = pd.read_csv("rmsdh_result/rmsdh_hingek_cnt_4.csv")
    small_fast_rmsdhk_2_postpro = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_2_postpro.csv"
    )
    small_fast_rmsdhk_3_postpro = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_3_postpro.csv"
    )
    small_fast_rmsdhk_4_postpro = pd.read_csv(
        "rmsdh_result/fast_rmsdh_hingek_cnt_4_postpro.csv"
    )

    fast_rmsdhk_2 = pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_2.csv")
    rmsdhk_2 = pd.read_csv("rmsdh_result/rmsdhk_more_data_2.csv")
    fast_rmsdhk_3 = pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_3.csv")
    rmsdhk_3 = pd.read_csv("rmsdh_result/rmsdhk_more_data_3.csv")
    fast_rmsdhk_4 = pd.read_csv("rmsdh_result/fast_rmsdhk_more_data_4.csv")
    rmsdhk_4 = pd.read_csv("rmsdh_result/rmsdhk_more_data_4.csv")
    fast_rmsdhk_postpro_2 = pd.read_csv(
        "rmsdh_result/fast_rmsdhk_more_data_2_pospro.csv"
    )
    fast_rmsdhk_postpro_3 = pd.read_csv(
        "rmsdh_result/fast_rmsdhk_more_data_3_pospro.csv"
    )
    fast_rmsdhk_postpro_4 = pd.read_csv(
        "rmsdh_result/fast_rmsdhk_more_data_4_pospro.csv"
    )
    more_monge_2 = preprocess_df(more_monge_2)
    rmsdhk_2 = preprocess_df(rmsdhk_2)
    rmsdhk_3 = preprocess_df(rmsdhk_3)
    rmsdhk_4 = preprocess_df(rmsdhk_4)
    fast_rmsdhk_2 = preprocess_df(fast_rmsdhk_2)
    fast_rmsdhk_3 = preprocess_df(fast_rmsdhk_3)
    fast_rmsdhk_4 = preprocess_df(fast_rmsdhk_4)
    fast_rmsdhk_postpro_2 = preprocess_df(fast_rmsdhk_postpro_2)
    fast_rmsdhk_postpro_3 = preprocess_df(fast_rmsdhk_postpro_3)
    fast_rmsdhk_postpro_4 = preprocess_df(fast_rmsdhk_postpro_4)
    for df in [
        rmsdhk_2,
        rmsdhk_3,
        rmsdhk_4,
        fast_rmsdhk_2,
        fast_rmsdhk_4,
        fast_rmsdhk_4,
        fast_rmsdhk_postpro_2,
        fast_rmsdhk_postpro_3,
        fast_rmsdhk_postpro_4,
    ]:
        print_exec_time(df)
    df_2 = calc_approximation_ratio_df(
        rmsdhk_2, fast_rmsdhk_2, more_monge_2=more_monge_2
    )
    df_3 = calc_approximation_ratio_df(
        rmsdhk_3, fast_rmsdhk_3, more_monge_2=more_monge_2
    )
    df_4 = calc_approximation_ratio_df(
        rmsdhk_4, fast_rmsdhk_4, more_monge_2=more_monge_2
    )
    df_2_postpro = calc_approximation_ratio_df(
        rmsdhk_2, fast_rmsdhk_postpro_2, more_monge_2=more_monge_2
    )
    df_3_postpro = calc_approximation_ratio_df(
        rmsdhk_3, fast_rmsdhk_postpro_3, more_monge_2=more_monge_2
    )
    df_4_postpro = calc_approximation_ratio_df(
        rmsdhk_4, fast_rmsdhk_postpro_4, more_monge_2=more_monge_2
    )
    small_df_2 = calc_approximation_ratio_df(
        small_rmsdhk_2, small_fast_rmsdhk_2, more_monge_2=small_monge_2
    )
    small_df_3 = calc_approximation_ratio_df(
        small_rmsdhk_3, small_fast_rmsdhk_3, more_monge_2=small_monge_2
    )
    small_df_4 = calc_approximation_ratio_df(
        small_rmsdhk_4, small_fast_rmsdhk_4, more_monge_2=small_monge_2
    )
    small_df_2_postpro = calc_approximation_ratio_df(
        small_rmsdhk_2, small_fast_rmsdhk_2_postpro, more_monge_2=small_monge_2
    )
    small_df_3_postpro = calc_approximation_ratio_df(
        small_rmsdhk_3, small_fast_rmsdhk_3_postpro, more_monge_2=small_monge_2
    )
    small_df_4_postpro = calc_approximation_ratio_df(
        small_rmsdhk_4, small_fast_rmsdhk_4_postpro, more_monge_2=small_monge_2
    )
    for df in [
        df_2,
        df_3,
        df_4,
        df_2_postpro,
        df_3_postpro,
        df_4_postpro,
        small_df_2,
        small_df_3,
        small_df_4,
        small_df_2_postpro,
        small_df_3_postpro,
        small_df_4_postpro,
    ]:
        print_approximation_ratio(df)
    plot_tmr(df_2, small_monge_2)
    calc_ans_shibuya(small_fast_rmsdhk_2, small_fast_rmsdhk_2_postpro)


def preprocess_df(df):
    ok_pdb = pd.read_csv("ok_pdb.csv")
    df = df[df["p_pdb_id"].isin(ok_pdb["p_pdb_id"])].reset_index()
    return df


def print_exec_time(df):
    print("Total computation time")
    print(df["exec_time (s)"].sum())
    print("Average computation time")
    print(df["exec_time (s)"].mean())


def calc_ans(detect_method_list, d=0):
    ans_list = []
    for i in range(len(expert_annotated_list)):
        exp = expert_annotated_list[i]
        detect = detect_method_list[i]
        true_hinge_indices = exp.split(" : ")
        detected_hinge_indices = detect.split(" : ")
        TP = 0
        FP = 0
        FN = 0
        if detected_hinge_indices != [""]:
            detected_ranges = [
                (int(label) - d, int(label) + d) for label in detected_hinge_indices
            ]
        else:
            detected_ranges = []
        for true in true_hinge_indices:
            if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
                TP += 1
            else:
                FN += 1
        for lower, upper in detected_ranges:
            if not any(lower <= int(true) <= upper for true in true_hinge_indices):
                FP += 1
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        ans_list.append(
            {"precision": precision, "recall": recall, "f_measure": f_measure}
        )
    return ans_list


def calc_ans_par(exp, detect, d=0):
    true_hinge_indices = exp.split(" : ")
    detected_hinge_indices = detect.split(" : ")
    TP = 0
    FP = 0
    FN = 0
    if detected_hinge_indices != [""]:
        detected_ranges = [
            (int(label) - d, int(label) + d) for label in detected_hinge_indices
        ]
    else:
        detected_ranges = []
    for true in true_hinge_indices:
        if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
            TP += 1
        else:
            FN += 1
    for lower, upper in detected_ranges:
        if not any(lower <= int(true) <= upper for true in true_hinge_indices):
            FP += 1
    return {"TP": TP, "FP": FP, "FN": FN}


def calc_approximation_ratio_df(df, heuristic_df, more_monge_2):
    df = pd.DataFrame(
        {
            "tmr": more_monge_2["monotonicity_rate_1"].to_list(),
            "RMSDhkD": (heuristic_df["RMSDh"] - df["RMSDh"]).to_list(),
            "approximation_ratio": (heuristic_df["RMSDh"] / df["RMSDh"]).to_list(),
            "RMSDhk": df["RMSDh"].to_list(),
            "heuristic_RMSDhk": heuristic_df["RMSDh"].to_list(),
        },
        index=more_monge_2.index,
    )
    return df


def print_approximation_ratio(df):
    print("Average approximation ratio")
    print(df["approximation_ratio"].mean())
    print("Max approximation ratio")
    print(df["approximation_ratio"].max())
    print("Min approximation ratio")
    print(df["approximation_ratio"].min())


def plot_tmr(df_2, small_monge_2):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 2))
    small_monge_2["monotonicity_rate_2"].hist(
        ax=axes[0],
        bins=[0.825, 0.85, 0.875, 0.900, 0.925, 0.950, 0.975, 1.00],
        edgecolor="none",
        color="Black",
    )
    axes[0].set_xlabel("Total monotonicity rate")
    axes[0].set_title("Distribution of TMR in the Shibuya 2008 dataset")
    axes[0].set_yticks([0, 1, 2, 3, 4])
    axes[0].set_xlim(0.825, 1.0)
    axes[0].set_ylim(0, 5)
    df_2["tmr"].hist(
        ax=axes[1],
        bins=[0.825, 0.85, 0.875, 0.900, 0.925, 0.950, 0.975, 1.00],
        edgecolor="none",
        color="Black",
    )
    axes[1].set_xlabel("Total monotonicity rate")
    axes[1].set_title("Distribution of TMR in the PAR 2020 dataset")
    axes[1].set_xlim(0.825, 1.0)
    axes[1].set_ylim(0, 40)
    plt.savefig("figures/distribution_of_tmr.svg", bbox_inches="tight", format="svg")
    plt.close()


def calc_ans_tp(detect_method_list, d=0):
    ans_list = []
    for i in range(len(expert_annotated_list)):
        exp = expert_annotated_list[i]
        detect = detect_method_list[i]
        true_hinge_indices = exp.split(" : ")
        detected_hinge_indices = detect.split(" : ")
        TP = 0
        FP = 0
        FN = 0
        if detected_hinge_indices != [""]:
            detected_ranges = [
                (int(label) - d, int(label) + d) for label in detected_hinge_indices
            ]
        else:
            detected_ranges = []
        for true in true_hinge_indices:
            if any(lower <= int(true) <= upper for (lower, upper) in detected_ranges):
                TP += 1
            else:
                FN += 1
        for lower, upper in detected_ranges:
            if not any(lower <= int(true) <= upper for true in true_hinge_indices):
                FP += 1
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        ans_list.append({"TP": TP, "FP": FP, "FN": FN})
    return ans_list


def calc_ans_shibuya(small_fast_rmsdhk_2, small_fast_rmsdhk_2_postpro, d=3):
    tmh_method_list = small_fast_rmsdhk_2["hinge_index"].tolist()
    tmhp_method_list = small_fast_rmsdhk_2_postpro["hinge_index"].tolist()
    method_f_measure_list = []
    for method, method_list in {
        "tmh": tmh_method_list,
        "tmhp": tmhp_method_list,
        "shibuya": shibuya_method_list,
        "flexprot": flexprot_list,
        "fatcat": fatcat_list,
        "nahal": nahal_list,
    }.items():
        TP = 0
        FP = 0
        FN = 0
        ans_df = pd.DataFrame(calc_ans_tp(method_list, d=d))
        TP = ans_df["TP"].sum()
        FP = ans_df["FP"].sum()
        FN = ans_df["FN"].sum()
        precision = TP / (TP + FP) if (TP + FP) != 0 else 0
        recall = TP / (TP + FN) if (TP + FN) != 0 else 0
        f_measure = (
            (2 * precision * recall) / (precision + recall)
            if (precision + recall) != 0
            else 0
        )
        method_f_measure_list.append(
            {
                "method": method,
                "d": d,
                "F-measure": f_measure,
                "Precision": precision,
                "Recall": recall,
            }
        )
    print(pd.DataFrame(method_f_measure_list))


if __name__ == "__main__":
    main()
