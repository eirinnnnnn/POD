#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


def metric_columns(names):
    pat = re.compile(r"^metric_([0-9]+)$")
    cols = []
    for name in names:
        m = pat.match(name)
        if m:
            cols.append((int(m.group(1)), name))
    return [name for _, name in sorted(cols)]


def read_dataset(path):
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    if not rows:
        raise SystemExit(f"empty dataset: {path}")
    metric_cols = metric_columns(rows[0].keys())
    M = len(metric_cols)
    metrics = np.array([[float(r[c]) for c in metric_cols] for r in rows], dtype=np.float64)
    labels = np.array([int(r["metric_argmin"]) for r in rows], dtype=np.int64)
    teacher_correct = np.array([int(r["teacher_correct"]) != 0 for r in rows], dtype=bool)
    correct = np.array([[int(r[f"correct_{i}"]) != 0 for i in range(M)] for r in rows], dtype=bool)
    return rows, metrics, labels, teacher_correct, correct


def read_predictions(path, M):
    pred_rows = list(csv.DictReader(open(path, encoding="utf-8")))
    if not pred_rows:
        raise SystemExit(f"empty prediction file: {path}")
    cols = [f"idx{i}" for i in range(M) if f"idx{i}" in pred_rows[0]]
    if not cols:
        raise SystemExit(f"prediction file has no idx columns: {path}")
    return np.array([[int(r[c]) for c in cols] for r in pred_rows], dtype=np.int64)


def topk_metrics(pred, metrics, labels, correct, basin, ks):
    out = []
    n, M = metrics.shape
    for k in ks:
        kk = min(k, pred.shape[1], M)
        top = pred[:, :kk]
        arg_hit = (top == labels[:, None]).any(axis=1).mean()
        basin_hit = basin[np.arange(n)[:, None], top].any(axis=1).mean()

        chosen = []
        for r in range(n):
            candidates = top[r]
            cand_metrics = metrics[r, candidates]
            chosen.append(candidates[int(np.argmin(cand_metrics))])
        chosen = np.asarray(chosen, dtype=np.int64)
        bler = 1.0 - correct[np.arange(n), chosen].mean()
        out.append((k, bler, basin_hit, arg_hit))
    return out


def describe_counter(counter, total, limit=10):
    return " ".join(f"{idx}:{cnt}({cnt/total:.3f})" for idx, cnt in counter.most_common(limit))


def bin_indices(basin_size, gaps):
    bins = []
    bins.append(("all", np.ones_like(basin_size, dtype=bool)))
    for name, mask in [
        ("basin_1", basin_size == 1),
        ("basin_2_4", (basin_size >= 2) & (basin_size <= 4)),
        ("basin_5_16", (basin_size >= 5) & (basin_size <= 16)),
        ("basin_17_63", (basin_size >= 17) & (basin_size <= 63)),
        ("basin_all", basin_size == basin_size.max()),
        ("gap_pos", gaps > 1e-12),
        ("gap_zero", gaps <= 1e-12),
    ]:
        if mask.any():
            bins.append((name, mask))
    return bins


def main():
    parser = argparse.ArgumentParser(description="Audit M-classifier datasets and prediction files.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--name", default="")
    parser.add_argument("--pred", action="append", nargs=2, metavar=("MODEL", "CSV"), default=[])
    parser.add_argument("--ks", default="1,2,4,8,16,32,64,128")
    parser.add_argument("--out", help="Optional report path.")
    args = parser.parse_args()

    rows, metrics, labels, teacher_correct, correct = read_dataset(args.dataset)
    n, M = metrics.shape
    ks = [int(x) for x in args.ks.split(",") if x.strip()]
    ks = [k for k in ks if k <= M]

    best = metrics.min(axis=1, keepdims=True)
    basin = np.abs(metrics - best) <= 1e-12
    basin_size = basin.sum(axis=1)
    sorted_metrics = np.sort(metrics, axis=1)
    gaps = sorted_metrics[:, 1] - sorted_metrics[:, 0] if M > 1 else np.zeros(n)

    lines = []
    title = args.name or Path(args.dataset).name
    lines.append(f"# {title}")
    lines.append(f"dataset,{args.dataset}")
    lines.append(f"samples,{n}")
    lines.append(f"M,{M}")
    lines.append(f"teacher_bler,{1.0 - teacher_correct.mean():.6f}")
    label_counter = Counter(labels.tolist())
    lines.append(f"unique_labels,{len(label_counter)}")
    lines.append(f"label_top10,{describe_counter(label_counter, n)}")
    lines.append(f"label_entropy,{-(sum((c/n)*math.log(max(c/n,1e-300)) for c in label_counter.values())):.6f}")
    lines.append(
        "basin_size,mean={:.6f},median={:.1f},p75={:.1f},p90={:.1f},max={},allM={}".format(
            float(basin_size.mean()),
            float(np.median(basin_size)),
            float(np.quantile(basin_size, 0.75)),
            float(np.quantile(basin_size, 0.90)),
            int(basin_size.max()),
            int((basin_size == M).sum()),
        )
    )
    lines.append(
        "gap,positive_frac={:.6f},mean={:.12g},median={:.12g},p90={:.12g},max={:.12g}".format(
            float((gaps > 1e-12).mean()),
            float(gaps.mean()),
            float(np.median(gaps)),
            float(np.quantile(gaps, 0.90)),
            float(gaps.max()),
        )
    )

    majority_order = [idx for idx, _ in label_counter.most_common()] + [i for i in range(M) if i not in label_counter]
    majority_pred = np.tile(np.array(majority_order, dtype=np.int64), (n, 1))
    pred_specs = [("majority", majority_pred)]
    for model_name, pred_path in args.pred:
        pred_specs.append((model_name, read_predictions(pred_path, M)))

    lines.append("")
    lines.append("## Prediction Summary")
    lines.append("model,unique_top1,top1_dist,top1_bler,top1_basin_hit,top1_argmin_hit")
    for model_name, pred in pred_specs:
        top1_counter = Counter(pred[:, 0].tolist())
        metric = topk_metrics(pred, metrics, labels, correct, basin, [1])[0]
        _, bler, basin_hit, arg_hit = metric
        lines.append(
            f"{model_name},{len(top1_counter)},{describe_counter(top1_counter, n, 8)},"
            f"{bler:.6f},{basin_hit:.6f},{arg_hit:.6f}"
        )

    lines.append("")
    lines.append("## Top-k")
    lines.append("model,k,bler,basin_hit,argmin_hit")
    for model_name, pred in pred_specs:
        for k, bler, basin_hit, arg_hit in topk_metrics(pred, metrics, labels, correct, basin, ks):
            lines.append(f"{model_name},{k},{bler:.6f},{basin_hit:.6f},{arg_hit:.6f}")

    lines.append("")
    lines.append("## Stratified Top1")
    lines.append("model,bin,count,teacher_bler,top1_bler,top1_basin_hit,top1_argmin_hit")
    for model_name, pred in pred_specs:
        top1 = pred[:, 0]
        for bin_name, mask in bin_indices(basin_size, gaps):
            count = int(mask.sum())
            if count == 0:
                continue
            top1_bler = 1.0 - correct[np.arange(n)[mask], top1[mask]].mean()
            top1_basin_hit = basin[np.arange(n)[mask], top1[mask]].mean()
            top1_arg_hit = (top1[mask] == labels[mask]).mean()
            teacher_bler = 1.0 - teacher_correct[mask].mean()
            lines.append(
                f"{model_name},{bin_name},{count},{teacher_bler:.6f},"
                f"{top1_bler:.6f},{top1_basin_hit:.6f},{top1_arg_hit:.6f}"
            )

    text = "\n".join(lines) + "\n"
    print(text)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
