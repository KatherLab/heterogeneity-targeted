# %%

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

rc_params = {
    "font.size": 10,  # 8
    "font.family": "DejaVu sans",
}

mpl.rcParams.update(rc_params)
mpl.style.use("bmh")


DataFrameColumn = str
"""The title of a Dataframe column"""

pos_class = "MUT"
neg_class = "WT"


def make_surv_df(
    df_clini: pd.DataFrame,
    duration: DataFrameColumn,
    event: DataFrameColumn,
    variate: DataFrameColumn,
) -> pd.DataFrame:
    """Return the processed dataframe for survival analysis"""
    df = df_clini[[duration, event, variate]].copy()
    # df_nas = df.isnull().sum()
    # print(f"df construction of {variate} | number of nulls {len(df_nas)}")
    df = df.dropna()
    return df


def coxph_entropy_pred(
    df_merge: pd.DataFrame, duration_col: DataFrameColumn, event: DataFrameColumn
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cph = CoxPHFitter()
    df_uni_entropy = make_surv_df(
        df_merge, duration_col, event, variate="entropy_binarized"
    )
    cph.fit(df_uni_entropy, duration_col=duration_col, event_col=event)
    df_results_univar = cph.summary

    # binarized patient-preds
    df_uni_pred = make_surv_df(df_merge, duration_col, event, "pred")
    cph = CoxPHFitter()
    cph.fit(df_uni_pred, duration_col=duration_col, event_col=event)
    df_results_univar = pd.concat([df_results_univar, cph.summary], axis=0)

    cph = CoxPHFitter()
    df_multi_entropy = df_merge[
        [duration_col, event, "entropy_binarized", "AGE", "GENDER", "STAGE"]
    ]
    df_multi_entropy.GENDER.replace({"male": 0, "female": 1}, inplace=True)
    df_multi_entropy = df_multi_entropy.dropna()
    cph.fit(df_multi_entropy, duration_col=duration_col, event_col=event)
    df_results_multivar = cph.summary

    cph = CoxPHFitter()
    df_multi_pred = df_merge[[duration_col, event, "pred", "AGE", "GENDER", "STAGE"]]
    df_multi_pred.GENDER.replace({"male": 0, "female": 1}, inplace=True)
    df_multi_pred = df_multi_entropy.dropna()
    cph.fit(df_multi_pred, duration_col=duration_col, event_col=event)
    df_results_multivar = pd.concat([df_results_multivar, cph.summary], axis=0)

    return df_uni_pred, df_results_univar, df_results_multivar


""" 
def univar_analysis(
    df_merge: pd.DataFrame, duration_col: DataFrameColumn, event: DataFrameColumn
) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    cph = CoxPHFitter()
    df_uni_entropy = make_surv_df(
        df_merge, duration_col, event, variate="entropy_binarized"
    )
    cph.fit(df_uni_entropy, duration_col=duration_col, event_col=event)
    df_results_univar = cph.summary

    # binarized patient-preds
    df_uni_pred = make_surv_df(df_merge, duration_col, event, "pred")
    cph = CoxPHFitter()
    cph.fit(df_uni_pred, duration_col=duration_col, event_col=event)
    df_results_univar = pd.concat([df_results_univar, cph.summary], axis=0)

    return df_uni_pred, df_results_univar


def multivar_analysis(
    df_merge: pd.DataFrame, duration_col: DataFrameColumn, event: DataFrameColumn
) -> pd.DataFrame:
 returns df uni pred and df results univar

    cph = CoxPHFitter()
    df_multi_entropy = df_merge[
        [duration_col, event, "entropy_binarized", "AGE", "GENDER", "STAGE"]
    ]
    df_multi_entropy.GENDER.replace({"male": 0, "female": 1}, inplace=True)
    df_multi_entropy = df_multi_entropy.dropna()
    cph.fit(df_multi_entropy, duration_col=duration_col, event_col=event)
    df_results_multivar = cph.summary

    cph = CoxPHFitter()
    df_multi_pred = df_merge[[duration_col, event, "pred", "AGE", "GENDER", "STAGE"]]
    df_multi_pred.GENDER.replace({"male": 0, "female": 1}, inplace=True)
    df_multi_pred = df_multi_entropy.dropna()
    cph.fit(df_multi_pred, duration_col=duration_col, event_col=event)
    df_results_multivar = pd.concat([df_results_multivar, cph.summary], axis=0)

    return df_results_multivar

"""


def coxph_combined_status(
    df_merge: pd.DataFrame, duration_col: DataFrameColumn, event_col: DataFrameColumn
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # deciphering the status
    df_merge = df_merge.copy()
    df_merge.loc[
        (df_merge["combined"] == 0) | (df_merge["combined"] == 2),
        f"{neg_class}_{pos_class}",
    ] = 0
    df_merge.loc[
        (df_merge["combined"] == 1) | (df_merge["combined"] == 3),
        f"{neg_class}_{pos_class}",
    ] = 1
    df_merge.loc[(df_merge["combined"] == 0) | (df_merge["combined"] == 1), "ITH"] = 0
    df_merge.loc[(df_merge["combined"] == 2) | (df_merge["combined"] == 3), "ITH"] = 1
    df_merge.loc[df_merge["combined"] == 0, f"{neg_class}_{pos_class}_LowITH"] = 0
    df_merge.loc[df_merge["combined"] == 1, f"{neg_class}_{pos_class}_LowITH"] = 1
    df_merge.loc[
        (df_merge["combined"] == 2) | (df_merge["combined"] == 3),
        f"{neg_class}_{pos_class}_LowITH",
    ] = np.nan
    df_merge.loc[df_merge["combined"] == 2, f"{neg_class}_{pos_class}_HighITH"] = 0
    df_merge.loc[df_merge["combined"] == 3, f"{neg_class}_{pos_class}_HighITH"] = 1
    df_merge.loc[
        (df_merge["combined"] == 0) | (df_merge["combined"] == 1),
        f"{neg_class}_{pos_class}_HighITH",
    ] = np.nan
    df_merge.loc[df_merge["combined"] == 0, f"Low_High_in{neg_class}"] = 0
    df_merge.loc[df_merge["combined"] == 2, f"Low_High_in{neg_class}"] = 1
    df_merge.loc[
        (df_merge["combined"] == 1) | (df_merge["combined"] == 3),
        f"Low_High_in{neg_class}",
    ] = np.nan
    df_merge.loc[df_merge["combined"] == 1, f"Low_High_in{pos_class}"] = 0
    df_merge.loc[df_merge["combined"] == 3, f"Low_High_in{pos_class}"] = 1
    df_merge.loc[
        (df_merge["combined"] == 0) | (df_merge["combined"] == 2),
        f"Low_High_in{pos_class}",
    ] = np.nan

    # univariate analysis
    cols = [
        f"{neg_class}_{pos_class}",
        "ITH",
        f"{neg_class}_{pos_class}_LowITH",
        f"{neg_class}_{pos_class}_HighITH",
        f"Low_High_in{neg_class}",
        f"Low_High_in{pos_class}",
    ]
    df_results_entropy_univar = pd.DataFrame()
    for col in cols:
        df_univar = make_surv_df(df_merge, duration_col, event_col, col)
        cph = CoxPHFitter()
        try:
            cph.fit(df_univar, duration_col=duration_col, event_col=event_col)
        except:
            continue
        df_results = cph.summary
        if df_results_entropy_univar.empty:
            df_results_entropy_univar = df_results.copy(deep=True)
        else:
            df_results_entropy_univar = pd.concat(
                [df_results_entropy_univar, df_results], axis=0
            )
    # multivariate analysis
    cols = [
        f"{neg_class}_{pos_class}",
        "ITH",
        f"{neg_class}_{pos_class}_LowITH",
        f"{neg_class}_{pos_class}_HighITH",
        f"Low_High_in{neg_class}",
        f"Low_High_in{pos_class}",
    ]
    df_results_entropy_multivar = pd.DataFrame()
    for col in cols:
        df_multivar = df_merge[[duration_col, event_col, col, "AGE", "GENDER", "STAGE"]]
        df_multivar.GENDER.replace({"male": 0, "female": 1}, inplace=True)
        # df_nas = df_multi.isnull().sum()
        df_multivar = df_multivar.dropna()
        cph = CoxPHFitter()
        try:
            cph.fit(df_multivar, duration_col=duration_col, event_col=event_col)
        except:
            continue
        df_results = cph.summary
        if df_results_entropy_multivar.empty:
            df_results_entropy_multivar = df_results.copy(deep=True)
        else:
            df_results_entropy_multivar = pd.concat(
                [df_results_entropy_multivar, df_results], axis=0
            )
    return df_results_entropy_univar, df_results_entropy_multivar


# %%


def run_survival_functions(
    target: str,
    path_entropy: str,
    path_clini_table: str,
    event_col: DataFrameColumn,
    duration_col: DataFrameColumn,
    path_survival: Path,
    path_svgs: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_entropy = pd.read_csv(path_entropy)
    df_clini = pd.read_csv(path_clini_table)

    df_merge = pd.merge(df_clini, df_entropy, on="PATIENT")
    # exclude stage 4 patients
    df_merge = df_merge[df_merge.STAGE != 4.0]
    n_stage_NaN = df_merge.STAGE.isnull().sum()
    df_merge = df_merge.dropna(subset=["STAGE"])

    # convert days to years & cap survival to 10 years
    df_merge["Fuyears"] = df_merge["Fudays"].map(lambda x: x / 360)
    df_merge.loc[df_merge["Fuyears"] >= 10, ["Fuyears", event_col]] = [10, 0]

    df_uni_pred, df_results_univar, df_results_multivar = coxph_entropy_pred(
        df_merge, duration_col, event_col
    )
    df_results_combined_univar, df_results_combined_multivar = coxph_combined_status(
        df_merge, duration_col, event_col
    )

    # %%

    # Kaplan Meier Curves for PATIENT PREDS
    kmf_pat = KaplanMeierFitter()
    ax = plt.subplot(111)
    T_low_pat = df_uni_pred[df_uni_pred[f"pred"] == 0][duration_col]
    E_low_pat = df_uni_pred[df_uni_pred[f"pred"] == 0][event_col]
    kmf_pat.fit(
        durations=T_low_pat,
        event_observed=E_low_pat,
        label=f"{target} Predicted WT (n={len(T_low_pat)})",
    )
    kmf_pat.plot_survival_function(ax=ax, color="navy", ci_show=False)  # C0
    T_high_pat = df_uni_pred[df_uni_pred[f"pred"] == 1][duration_col]
    E_high_pat = df_uni_pred[df_uni_pred[f"pred"] == 1][event_col]
    kmf_pat.fit(
        durations=T_high_pat,
        event_observed=E_high_pat,
        label=f"{target} Predicted MUT (n={len(T_high_pat)})",
    )
    kmf_pat.plot_survival_function(ax=ax, color="sienna", ci_show=False)
    km_results_patientPreds = logrank_test(
        T_low_pat, T_high_pat, E_low_pat, E_high_pat, alpha=0.99
    )

    ax.set_xlabel("Years", fontsize=15)
    ax.set_ylabel("Estimated Probability of Survival", fontsize=15)
    ax.set_title(
        f"log-rank p={km_results_patientPreds.p_value:.2e}", fontsize=15, y=1.02
    )
    ax.set_facecolor("#f9f9f9")
    plt.savefig(
        str(path_survival / f"KM_patient_preds_{target}_{event_col}.png"),
        dpi=1000,
        bbox_inches="tight",
    )
    plt.savefig(
        str(path_svgs / f"KM_patient_preds_{target}_{event_col}.svg"),
        format="svg",
        dpi=1000,
        bbox_inches="tight",
    )
    plt.clf()
    ax.cla()
    # %%
    df_four = df_merge[[duration_col, event_col, "combined"]]
    df_four = df_four.dropna()
    kmf_four = KaplanMeierFitter()
    ax = plt.subplot(111)
    # ITH Low, Predicted Low
    if not df_four[df_four.combined == 0].empty:
        T_lowE_lowP = df_four[df_four.combined == 0][duration_col]
        E_lowE_lowP = df_four[df_four.combined == 0][event_col]
        kmf_four.fit(
            durations=T_lowE_lowP,
            event_observed=E_lowE_lowP,
            label=f"{target} ITH Low & Predicted {neg_class} ({len(T_lowE_lowP)})",
        )
        kmf_four.plot_survival_function(ax=ax, color="navy", ci_show=False)  # 'C9'
    # ITH low, predicted high
    if not df_four[df_four.combined == 1].empty:
        T_lowE_highP = df_four[df_four.combined == 1][duration_col]
        E_lowE_highP = df_four[df_four.combined == 1][event_col]
        kmf_four.fit(
            durations=T_lowE_highP,
            event_observed=E_lowE_highP,
            label=f"{target} ITH Low & Predicted {pos_class} ({len(T_lowE_highP)})",
        )
        kmf_four.plot_survival_function(ax=ax, color="C5", ci_show=False)
    # ITH high, predicted low
    if not df_four[df_four.combined == 2].empty:
        T_highE_lowP = df_four[df_four.combined == 2][duration_col]
        E_highE_lowP = df_four[df_four.combined == 2][event_col]
        kmf_four.fit(
            durations=T_highE_lowP,
            event_observed=E_highE_lowP,
            label=f"{target} ITH High & Predicted {neg_class} ({len(T_highE_lowP)})",
        )
        kmf_four.plot_survival_function(ax=ax, color="sienna", ci_show=False)
    # ITH high, predicted High
    if not df_four[df_four.combined == 3].empty:
        T_highE_highP = df_four[df_four.combined == 3][duration_col]
        E_highE_highP = df_four[df_four.combined == 3][event_col]
        kmf_four.fit(
            durations=T_highE_highP,
            event_observed=E_highE_highP,
            label=f"{target} ITH High & Predicted {pos_class} ({len(T_highE_highP)})",
        )
        kmf_four.plot_survival_function(ax=ax, color="darkgreen", ci_show=False)
    # %%
    result_multi = multivariate_logrank_test(
        df_four[duration_col], df_four["combined"], df_four[event_col]
    )
    ax.set_facecolor("#f9f9f9")
    ax.set_xlabel("Years", fontsize=15)
    ax.set_ylabel("Estimated Probability of Survival", fontsize=15)
    ax.set_title(f"log-rank p={result_multi.p_value:.2e}", fontsize=15, y=1.02)

    plt.savefig(
        str(path_survival / f"KM_combined_4_{target}_{event_col}.png"),
        dpi=1000,
        bbox_inches="tight",
    )
    plt.savefig(
        str(path_svgs / f"KM_combined_4_{target}_{event_col}.svg"),
        format="svg",
        dpi=1000,
        bbox_inches="tight",
    )
    plt.clf()
    ax.cla()

    df_results_univar.rename(
        columns={"entropy_binarized": f"{target}_ITH"}, inplace=True
    )
    df_results_multivar.rename(
        columns={"entropy_binarized": f"{target}_ITH"}, inplace=True
    )

    # %%
    print(f"Number of patients with NaN stage values: {n_stage_NaN}")

    return (
        df_results_univar,
        df_results_multivar,
        df_results_combined_univar,
        df_results_combined_multivar,
    )
