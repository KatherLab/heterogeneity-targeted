import json

import numpy as np
import pandas as pd


def comp_entropy(path_tile_preds: str, path_patient_pred: str) -> pd.DataFrame:
    df_preds = pd.read_csv(path_tile_preds)
    df_pooled_preds = pd.read_csv(path_patient_pred)

    cls_preds1, cls_preds2 = list(df_preds.columns)[1:3]
    cls_name = cls_preds1.split("_")[1]
    df_preds[f"is{cls_name}"] = df_preds[cls_preds1] > df_preds[cls_preds2]

    patient_names = df_preds["PATIENT"].unique()
    data = {"PATIENT": [], "entropy": []}
    for patient in patient_names:
        df_patient_preds = df_preds[df_preds["PATIENT"] == patient]
        len_df_pat = len(df_patient_preds)

        dict_is_high = json.loads(
            df_patient_preds[f"is{cls_name}"].value_counts().to_json()
        )
        sum = 0
        for key in dict_is_high:
            prob = dict_is_high[key] / len_df_pat
            sum += prob * np.log2(prob)
        s_ent = -sum

        data["PATIENT"].append(patient)
        data["entropy"].append(s_ent)

    df_entropy = pd.DataFrame.from_dict(data)

    df_merge = pd.merge(df_entropy, df_pooled_preds, on="PATIENT")
    median_val = df_merge["entropy"].median()

    df_merge["entropy_binarized"] = df_merge["entropy"].map(
        lambda x: 0 if x < median_val else 1
    )
    df_merge.loc[
        df_merge["entropy_binarized"] == 0, "entropy_binarized_decipher"
    ] = "Low ITH"
    df_merge.loc[
        df_merge["entropy_binarized"] == 1, "entropy_binarized_decipher"
    ] = "High ITH"

    # %%
    df_merge["pred_decipher"] = df_merge["pred"].copy(deep=True)
    df_merge[f"pred"].replace({"MUT": 1, "WT": 0}, inplace=True)

    df_merge["combined"] = df_merge["entropy_binarized"].astype(str) + df_merge[
        f"pred"
    ].astype(str)

    # Map the combinations to integers
    mapping = {"00": 0, "01": 1, "10": 2, "11": 3}

    df_merge["combined"] = df_merge["combined"].map(mapping)
    df_merge.loc[
        df_merge["combined"] == 0, "combined_decipher"
    ] = "ITH Low & Predicted WT"
    df_merge.loc[
        df_merge["combined"] == 1, "combined_decipher"
    ] = "ITH Low & Predicted MUT"
    df_merge.loc[
        df_merge["combined"] == 2, "combined_decipher"
    ] = "ITH High & Predicted WT"
    df_merge.loc[
        df_merge["combined"] == 3, "combined_decipher"
    ] = "ITH High & Predicted MUT"

    df_entropy = df_merge[
        [
            "PATIENT",
            "entropy",
            "entropy_binarized",
            "entropy_binarized_decipher",
            "pred",
            "pred_decipher",
            "combined",
            "combined_decipher",
        ]
    ]

    return df_entropy
