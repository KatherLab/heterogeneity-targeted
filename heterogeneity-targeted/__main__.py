#!/usr/bin/env python3
import argparse
from pathlib import Path

from .measure_entropy import comp_entropy
from .survival_analyses import run_survival_functions


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Visualize the t-SNE representation of images"
    )

    parser.add_argument(
        "--target",
        type=str,
        required=True,
        help="The target genetic alteration to be analyzed",
    )
    parser.add_argument(
        "--path-tile-preds",
        type=str,
        required=True,
        help="The path to tile-level predictions",
    )
    parser.add_argument(
        "--path-patient-pred",
        type=str,
        required=True,
        help="The path to patient-level predictions",
    )
    parser.add_argument(
        "--path-clini-table", type=str, required=True, help="The path to clinical table"
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        type=str,
        required=False,
        default="out",
        help="Path to folder where the results will be saved",
    )
    parser.add_argument(
        "--duration-col",
        type=str,
        required=False,
        help="The column name indicating survival time in days",
    )
    parser.add_argument(
        "--event-col",
        type=str,
        required=False,
        help="The column name indicating the event",
    )
    args = parser.parse_args()

    path_entropy_scores = Path(args.output_directory) / "scores_heterogeneity"
    path_entropy_scores.mkdir(exist_ok=True, parents=True)
    path_entropy_csv = path_entropy_scores / f"{args.target}__patient_to_entropy.csv"
    if not path_entropy_scores.exists():
        df_entropy = comp_entropy(args.path_tile_preds, args.path_patient_pred)
        df_entropy.to_csv(path_entropy_csv)

    path_survival = (
        Path(args.output_directory) / "Results" / args.target / args.event_col
    )
    path_survival.mkdir(exist_ok=True, parents=True)
    path_svgs = path_survival / "svg_images"
    path_svgs.mkdir(exist_ok=True, parents=True)

    (
        df_results_univar,
        df_results_multivar,
        df_results_combined_univar,
        df_results_combined_multivar,
    ) = run_survival_functions(
        args.target,
        path_entropy_csv,
        args.path_clini_table,
        args.event_col,
        args.duration_col,
        path_survival,
        path_svgs,
    )

    df_results_univar.to_csv(
        path_survival / f"{args.target}_CoxPH_univariate_{args.event_col}.csv"
    )
    df_results_multivar.to_csv(
        path_survival / f"{args.target}_CoxPH_multivariate_{args.event_col}.csv"
    )
    df_results_combined_univar.to_csv(
        path_survival / f"{args.target}_CoxPH_COMBINED_univariate_{args.event_col}.csv"
    )
    df_results_combined_multivar.to_csv(
        path_survival
        / f"{args.target}_CoxPH_COMBINED_multivariate_{args.event_col}.csv"
    )


if __name__ == "__main__":
    main()
