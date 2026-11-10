"""Interactive multi-EIC fast LC-MS analysis used for metallaphotoredox sulfonamidation.

Select one exported CSV/TXT file per EIC. The program aligns
cycles, finds each EIC apex in its assigned retention-time window, and exports as CSV tables. Requires numpy.
"""

import csv
import math
import re
from pathlib import Path
from tkinter import Tk, filedialog

import numpy as np


def ask_text(question, default):
    answer = input(f"{question} [{default}]: ").strip()
    return answer if answer else str(default)


def ask_float(question, default):
    while True:
        try:
            return float(ask_text(question, default))
        except ValueError:
            print("Please enter a number.")


def ask_integer(question, default):
    while True:
        try:
            return int(ask_text(question, default))
        except ValueError:
            print("Please enter a whole number.")


def choose_files():
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    filenames = filedialog.askopenfilenames(
        title="Choose one exported EIC file for each analyte",
        filetypes=[
            ("Chromatographic exports", "*.csv *.txt *.tsv"),
            ("All files", "*.*"),
        ],
    )

    root.destroy()

    if not filenames:
        raise SystemExit("No EIC files were selected.")

    return [Path(item) for item in filenames]


def detect_delimiter(path):
    preview = path.read_text(errors="replace")[:5000]
    return "\t" if preview.count("\t") >= preview.count(",") else ","


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_") or "analyte"


def quadratic_apex(time_min, signal, index):
    """Three-point quadratic refinement of the measured maximum."""
    if index == 0 or index == len(time_min) - 1:
        return float(time_min[index]), float(signal[index])

    x = time_min[index - 1:index + 2]
    y = signal[index - 1:index + 2]
    a, b, c = np.polyfit(x, y, 2)

    if a >= 0:
        return float(time_min[index]), float(signal[index])

    peak_x = -b / (2 * a)

    if not x[0] <= peak_x <= x[-1]:
        return float(time_min[index]), float(signal[index])

    return float(peak_x), float(a * peak_x**2 + b * peak_x + c)


def write_csv(path, headers, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def load_trace(
    path,
    header_rows,
    time_column,
    signal_column,
    start_min,
    end_min,
):
    raw = np.genfromtxt(
        path,
        delimiter=detect_delimiter(path),
        skip_header=header_rows,
        usecols=(time_column, signal_column),
        invalid_raise=False,
    )

    raw = np.atleast_2d(raw)
    raw = raw[np.isfinite(raw[:, 0]) & np.isfinite(raw[:, 1])]

    selected = (raw[:, 0] >= start_min) & (raw[:, 0] <= end_min)

    if not np.any(selected):
        raise ValueError(
            f"No numeric data in {path.name} within "
            f"{start_min}–{end_min} min."
        )

    return raw[selected, 0] - start_min, raw[selected, 1]


def main():
    print("\nFast-LC multi-EIC analysis\n")

    paths = choose_files()

    print("Selected EIC files:")
    for path in paths:
        print(f"  {path.name}")

    print("\nUse the same settings for all selected files.")

    header_rows = ask_integer("Header rows to skip", 1)
    time_column = ask_integer(
        "Time column (zero-based; first column = 0)", 0
    )
    signal_column = ask_integer(
        "Signal column (zero-based; second column = 1)", 1
    )
    start_min = ask_float("Start of analysis window, min", 0.0)
    end_min = ask_float("End of analysis window, min", 30.0)
    cycle_s = ask_float("Cycle period, s", 40.0)

    default_cycles = round((end_min - start_min) * 60 / cycle_s)
    n_cycles = ask_integer("Number of injection cycles", default_cycles)

    display_step_min = ask_float(
        "Displayed run-time interval, min", 1.0
    )

    print(
        "\nFor each EIC, enter its displayed analyte name and the "
        "expected RT window within one cycle."
    )

    eics = []

    for path in paths:
        name = ask_text(f"Analyte name for {path.name}", path.stem)

        left_s = ask_float(
            f"{name}: left RT-window edge, s", 0.0
        )
        right_s = ask_float(
            f"{name}: right RT-window edge, s",
            math.floor(cycle_s * 10) / 10,
        )

        if right_s <= left_s:
            raise SystemExit(f"Invalid RT window for {name}.")

        time_min, signal = load_trace(
            path,
            header_rows,
            time_column,
            signal_column,
            start_min,
            end_min,
        )

        eics.append(
            {
                "name": name,
                "path": path,
                "left_s": left_s,
                "right_s": right_s,
                "time_min": time_min,
                "signal": signal,
            }
        )

    for eic in eics:
        rt_values = []
        heights = []

        for cycle in range(n_cycles):
            cycle_start = cycle * cycle_s / 60

            inside = (
                (
                    eic["time_min"]
                    >= cycle_start + eic["left_s"] / 60
                )
                & (
                    eic["time_min"]
                    <= cycle_start + eic["right_s"] / 60
                )
            )

            indices = np.flatnonzero(inside)

            if len(indices) == 0:
                raise SystemExit(
                    f"No data for {eic['name']} in injection {cycle + 1}."
                )

            apex_index = indices[np.argmax(eic["signal"][indices])]

            apex_time, apex_height = quadratic_apex(
                eic["time_min"],
                eic["signal"],
                int(apex_index),
            )

            rt_values.append(apex_time * 60 - cycle * cycle_s)
            heights.append(apex_height)

        eic["rt"] = np.asarray(rt_values)
        eic["height"] = np.asarray(heights)

    output_dir = paths[0].parent / "fast_lc_eic_analysis_output"
    output_dir.mkdir(exist_ok=True)

    names = [eic["name"] for eic in eics]
    run_times = np.arange(n_cycles) * display_step_min

    write_csv(
        output_dir / "peak_heights.csv",
        ["normalized_run_time_min", *names],
        [
            [run_times[i], *[eic["height"][i] for eic in eics]]
            for i in range(n_cycles)
        ],
    )

    write_csv(
        output_dir / "retention_times.csv",
        ["normalized_run_time_min", *names],
        [
            [run_times[i], *[eic["rt"][i] for eic in eics]]
            for i in range(n_cycles)
        ],
    )

    rt_summary = []
    height_summary = []

    for eic in eics:
        rt = eic["rt"]
        height = eic["height"]

        rt_summary.append(
            [
                eic["name"],
                np.mean(rt),
                np.std(rt, ddof=1),
                100 * np.std(rt, ddof=1) / np.mean(rt),
            ]
        )

        height_summary.append(
            [
                eic["name"],
                np.mean(height),
                np.std(height, ddof=1),
                100 * np.std(height, ddof=1) / np.mean(height),
            ]
        )

    write_csv(
        output_dir / "retention_time_summary.csv",
        ["analyte", "mean_rt_s", "sd_rt_s", "rt_rsd_percent"],
        rt_summary,
    )

    write_csv(
        output_dir / "peak_height_summary.csv",
        [
            "analyte",
            "mean_peak_height",
            "sd_peak_height",
            "peak_height_rsd_percent",
        ],
        height_summary,
    )

    raw_rows = []

    for eic in eics:
        raw_rows.extend(
            [
                [time, eic["name"], signal]
                for time, signal in zip(
                    eic["time_min"],
                    eic["signal"],
                )
            ]
        )

    write_csv(
        output_dir / "normalized_raw_eic_traces.csv",
        ["normalized_run_time_min", "analyte", "raw_signal"],
        raw_rows,
    )

    rt_grid = np.round(
        np.arange(
            0,
            math.floor(cycle_s * 10) / 10 + 0.01,
            0.1,
        ),
        10,
    )

    run_grid = np.round(
        np.arange(
            0,
            (n_cycles - 1) * display_step_min + 0.01,
            display_step_min,
        ),
        10,
    )

    for eic in eics:
        matrix = [
            [
                rt_s,
                *[
                    np.interp(
                        run_min + rt_s / 60,
                        eic["time_min"],
                        eic["signal"],
                    )
                    for run_min in run_grid
                ],
            ]
            for rt_s in rt_grid
        ]

        write_csv(
            output_dir / f"heatmap_{safe_name(eic['name'])}.csv",
            ["retention_time_s", *run_grid],
            matrix,
        )

    with open(
        output_dir / "analysis_settings.txt",
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"Analysis window: {start_min}–{end_min} min\n"
        )
        handle.write(f"Cycle period: {cycle_s} s\n")
        handle.write(f"Cycles: {n_cycles}\n")
        handle.write(
            "Processing: maximum within analyte-specific RT window; "
            "three-point quadratic apex refinement; linear "
            "interpolation for heatmap display only.\n"
        )
        handle.write(
            "No baseline subtraction, smoothing, or intensity "
            "normalization applied.\n\n"
        )
        handle.write("EIC files and windows:\n")

        for eic in eics:
            handle.write(
                f"{eic['name']}: {eic['path']} | "
                f"{eic['left_s']}–{eic['right_s']} s\n"
            )

    print(f"\nAnalysis complete. Files were saved to:\n{output_dir}\n")

    print("Retention-time summary:")
    for row in rt_summary:
        print(
            f"  {row[0]}: {row[1]:.3f} ± {row[2]:.3f} s "
            f"({row[3]:.2f}% RSD)"
        )

    print("\nPeak-height summary:")
    for row in height_summary:
        print(
            f"  {row[0]}: {row[1]:.2f} ± {row[2]:.2f} "
            f"({row[3]:.2f}% RSD)"
        )


if __name__ == "__main__":
    main()