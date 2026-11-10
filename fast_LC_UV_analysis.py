"""Fast LC-UV injection repeatability analysis used for thiol-ene.

A file-picker opens first. The defaults reproduce the
thiol–ene repeatability analysis. Press Enter at each prompt to accept them.

Requires: Python 3 and numpy (install in terminal with: py -m pip install numpy)
"""

import csv
import math
from pathlib import Path
from tkinter import Tk, filedialog

import numpy as np


DEFAULT_PEAKS = [
    ("p-Toluidine", 0.60, 1.25),
    ("Product", 1.25, 1.85),
    ("4-Chlorostyrene", 3.45, 4.55),
]


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


def choose_file():
    root = Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    filename = filedialog.askopenfilename(
        title="Choose exported chromatographic time-series data",
        filetypes=[
            ("Chromatographic exports", "*.csv *.txt *.tsv"),
            ("All files", "*.*"),
        ],
    )
    root.destroy()

    if not filename:
        raise SystemExit("No data file was selected.")

    return Path(filename)


def detect_delimiter(path):
    preview = path.read_text(errors="replace")[:5000]
    return "\t" if preview.count("\t") >= preview.count(",") else ","


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


def main():
    print("\nFast-LC cycle analysis\n")

    source = choose_file()
    print(f"Selected file: {source.name}\n")

    header_rows = ask_integer("Header rows to skip", 2)
    time_column = ask_integer(
        "Time column (zero-based; first column = 0)", 0
    )
    signal_column = ask_integer(
        "Signal column (zero-based; second column = 1)", 1
    )

    start_min = ask_float("Start of analysis window, min", 0.0)
    end_min = ask_float("End of analysis window, min", 3.0)
    cycle_s = ask_float("Cycle period, s", 6.017234482758621)

    default_cycles = round((end_min - start_min) * 60 / cycle_s)
    n_cycles = ask_integer("Number of injection cycles", default_cycles)
    display_step_min = ask_float(
        "Displayed run-time interval, min", 0.1
    )

    print(
        "\nPeak windows are retention times within one injection cycle "
        "(seconds)."
    )

    peaks = []
    for default_name, default_left, default_right in DEFAULT_PEAKS:
        name = ask_text("Peak name", default_name)
        left = ask_float(f"{name}: left edge, s", default_left)
        right = ask_float(f"{name}: right edge, s", default_right)

        if right <= left:
            raise SystemExit(
                f"Invalid window for {name}: right edge must exceed left edge."
            )

        peaks.append((name, left, right))

    delimiter = detect_delimiter(source)

    raw = np.genfromtxt(
        source,
        delimiter=delimiter,
        skip_header=header_rows,
        usecols=(time_column, signal_column),
        invalid_raise=False,
    )
    raw = np.atleast_2d(raw)

    valid = np.isfinite(raw[:, 0]) & np.isfinite(raw[:, 1])
    raw = raw[valid]

    time_min = raw[:, 0]
    signal = raw[:, 1]

    selected = (time_min >= start_min) & (time_min <= end_min)
    time_min = time_min[selected] - start_min
    signal = signal[selected]

    if len(time_min) == 0:
        raise SystemExit(
            "No numeric data were found in the selected analysis window."
        )

    peak_data = {}

    for name, left_s, right_s in peaks:
        rt_values = []
        height_values = []

        for cycle in range(n_cycles):
            cycle_start = cycle * cycle_s / 60

            inside = (
                (time_min >= cycle_start + left_s / 60)
                & (time_min <= cycle_start + right_s / 60)
            )
            indices = np.flatnonzero(inside)

            if len(indices) == 0:
                raise SystemExit(
                    f"No data found for {name} in injection {cycle + 1}."
                )

            apex_index = indices[np.argmax(signal[indices])]
            apex_time, apex_height = quadratic_apex(
                time_min, signal, int(apex_index)
            )

            rt_values.append(apex_time * 60 - cycle * cycle_s)
            height_values.append(apex_height)

        peak_data[name] = {
            "rt": np.asarray(rt_values),
            "height": np.asarray(height_values),
        }

    output_dir = source.parent / "fast_lc_analysis_output"
    output_dir.mkdir(exist_ok=True)

    names = [item[0] for item in peaks]
    run_times = np.arange(n_cycles) * display_step_min

    write_csv(
        output_dir / "peak_heights.csv",
        ["normalized_run_time_min", *names],
        [
            [run_times[i], *[peak_data[name]["height"][i] for name in names]]
            for i in range(n_cycles)
        ],
    )

    write_csv(
        output_dir / "retention_times.csv",
        ["normalized_run_time_min", *names],
        [
            [run_times[i], *[peak_data[name]["rt"][i] for name in names]]
            for i in range(n_cycles)
        ],
    )

    rt_summary = []
    height_summary = []

    for name in names:
        rt = peak_data[name]["rt"]
        height = peak_data[name]["height"]

        rt_summary.append(
            [
                name,
                np.mean(rt),
                np.std(rt, ddof=1),
                100 * np.std(rt, ddof=1) / np.mean(rt),
            ]
        )

        height_summary.append(
            [
                name,
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

    write_csv(
        output_dir / "normalized_raw_trace.csv",
        ["normalized_run_time_min", "raw_signal"],
        zip(time_min, signal),
    )

    rt_grid = np.round(
        np.arange(0, math.floor(cycle_s * 10) / 10 + 0.01, 0.1), 10
    )
    run_grid = np.round(
        np.arange(
            0,
            (n_cycles - 1) * display_step_min + 0.01,
            display_step_min,
        ),
        10,
    )

    heatmap_rows = [
        [
            rt_s,
            *[
                np.interp(run_min + rt_s / 60, time_min, signal)
                for run_min in run_grid
            ],
        ]
        for rt_s in rt_grid
    ]

    write_csv(
        output_dir / "heatmap_matrix.csv",
        ["retention_time_s", *run_grid],
        heatmap_rows,
    )

    with open(
        output_dir / "analysis_settings.txt", "w", encoding="utf-8"
    ) as handle:
        handle.write(f"Source file: {source}\n")
        handle.write(f"Analysis window: {start_min}–{end_min} min\n")
        handle.write(f"Cycle period: {cycle_s} s\n")
        handle.write(f"Cycles: {n_cycles}\n")
        handle.write("Peak windows (s):\n")

        for name, left_s, right_s in peaks:
            handle.write(f"  {name}: {left_s}–{right_s}\n")

        handle.write(
            "Processing: maximum within window; three-point quadratic "
            "apex refinement; linear interpolation for heatmap only.\n"
        )
        handle.write(
            "No baseline subtraction, smoothing, or intensity "
            "normalization applied.\n"
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