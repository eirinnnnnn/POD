#!/usr/bin/env python3
"""model_out variant: same SET-LEVEL MLP regressor as
adaptive_path_selector/pm/ebch_m7t10/train_adaptive_m.py (predicts
m_required_basin -- the smallest pruning width guaranteed to include a
correctly-decoding branch), but the input is the trace_learnt_path_selector
model's OWN dim-M output (a sorted per-branch probability/score vector)
instead of pm's raw sorted checkpoint pm_min vector.

Input/target come straight from mclass_adaptive_m_dataset_model_out's CSV:
    sample, M, m_required, full_ped_correct, m_required_basin,
    any_basin_exists, t{T}_trace_prob_sorted_1..M (or
    t{T}_trace_score_sorted_1..M for a direction_largest=false, i.e.
    log_gap-trained, trace model)

Run with --no-log-input: pm's log1p(pm_min) transform exists because raw
pm_min is an unbounded, heavy-tailed linear metric; trace_prob is already
a bounded [0,1] sigmoid output, so log-compressing it is the wrong
transform here, not a neutral no-op.

Self-contained: trains AND exports a plain-text weight file in one step,
same format as the pm sibling.
"""
import argparse
from pathlib import Path

import numpy as np


def load_dataset(path):
    with open(path, newline="", encoding="utf-8") as f:
        names = f.readline().strip().split(",")
    col = {name: i for i, name in enumerate(names)}
    raw = np.loadtxt(path, delimiter=",", skiprows=1, dtype=np.float64)
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)

    M = int(raw[0, col["M"]])
    feature_cols = [n for n in names if n.startswith("t") and "_sorted_" in n]
    feature_cols.sort(key=lambda n: int(n.rsplit("_", 1)[1]))
    if len(feature_cols) != M:
        raise ValueError(f"expected {M} sorted-feature columns, found {len(feature_cols)}")
    feature_idx = [col[n] for n in feature_cols]

    X = raw[:, feature_idx].astype(np.float32)
    # m_required_basin: smallest m such that top-m contains ANY correctly-
    # decoding branch (not necessarily full_ped's specific metric-argmin
    # choice) -- confirmed to reproduce mclass_bler_sim's real static_bler(k)
    # essentially exactly at every k, unlike the stricter single-teacher
    # m_required (see NOTES.md).
    y = raw[:, col["m_required_basin"]].astype(np.float32)
    full_ped_correct = raw[:, col["any_basin_exists"]].astype(np.float32)
    return X, y, full_ped_correct, M, feature_cols


def log_transform_inputs(X):
    # Not the natural transform for a trace_prob input (already bounded
    # [0,1]) -- kept only so --no-log-input has a symmetric opposite to
    # disable; pm's version applies this by default (raw pm_min is
    # unbounded/heavy-tailed), model_out should always be run WITHOUT it.
    return np.log1p(X)


def standardize(train, val):
    mean = train.mean(axis=0, keepdims=True)
    std = train.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return (train - mean) / std, (val - mean) / std, mean.reshape(-1), std.reshape(-1)


def write_vec(f, name, arr):
    arr = np.asarray(arr, dtype=np.float64).reshape(-1)
    f.write(f"{name} {len(arr)}\n")
    f.write(" ".join(repr(float(x)) for x in arr) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True)
    parser.add_argument("--out", required=True, help="plain-text weight file path")
    parser.add_argument("--save-pt", default=None, help="optional torch checkpoint path (for resuming/inspection)")
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--under-penalty", type=float, default=1.0,
                         help="loss weight multiplier when prediction < target (under-predicting m "
                              "risks missing full_ped's answer; over-predicting only costs compute). "
                              "1.0 = plain symmetric MSE.")
    parser.add_argument("--no-log-input", action="store_true",
                         help="skip the log1p(pm_min) input transform (on by default)")
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, y, full_ped_correct, M, feature_cols = load_dataset(args.data)
    log_input = not args.no_log_input
    if log_input:
        X = log_transform_inputs(X)
    n = X.shape[0]
    idx = np.arange(n)
    rng.shuffle(idx)
    val_count = max(1, int(n * args.val_frac))
    val_idx, train_idx = idx[:val_count], idx[val_count:]

    X_train, X_val, x_mean, x_std = standardize(X[train_idx], X[val_idx])
    y_mean = float(y[train_idx].mean())
    y_std = float(y[train_idx].std())
    if y_std < 1e-6:
        y_std = 1.0
    y_train = (y[train_idx] - y_mean) / y_std
    y_val = (y[val_idx] - y_mean) / y_std

    model = nn.Sequential(
        nn.Linear(M, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    loader = DataLoader(
        TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train)),
        batch_size=args.batch_size, shuffle=True,
    )
    Xv = torch.from_numpy(X_val)
    yv_raw = torch.from_numpy(y[val_idx])

    def asym_mse(pred, target):
        err = pred - target
        w = torch.where(err < 0, args.under_penalty, 1.0)
        return (w * err * err).mean()

    import copy
    best_mae = float("inf")
    best_epoch = 0
    best_state = None

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss, total = 0.0, 0
        for xb, yb in loader:
            opt.zero_grad()
            pred = model(xb).reshape(-1)
            loss = asym_mse(pred, yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * xb.shape[0]
            total += xb.shape[0]

        model.eval()
        with torch.no_grad():
            pred_v_raw = model(Xv).reshape(-1) * y_std + y_mean
            mae = float((pred_v_raw - yv_raw).abs().mean().item())
            under = float(((pred_v_raw - yv_raw) < 0).float().mean().item())
            # what would actually happen if we used ceil(pred) as the pruning
            # width: miss iff ceil(pred) < true m_required (and full_ped
            # itself was going to succeed -- if it wasn't, no m saves it)
            fpc_v = torch.from_numpy(full_ped_correct[val_idx])
            used_m = torch.clamp(torch.ceil(pred_v_raw), 1, M)
            miss = float(((used_m < yv_raw) & (fpc_v > 0.5)).float().mean().item())
            mean_m_used = float(used_m.mean().item())
        print(f"epoch={epoch} train_loss={total_loss/total:.6f} val_mae={mae:.4f} "
              f"val_under_rate={under:.4f} val_miss_rate={miss:.4f} val_mean_m={mean_m_used:.2f} (M={M})")

        if mae < best_mae:
            best_mae = mae
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())

    print(f"selecting best epoch={best_epoch} (val_mae={best_mae:.4f}) for export")
    model.load_state_dict(best_state)

    if args.save_pt:
        import torch as _torch
        Path(args.save_pt).parent.mkdir(parents=True, exist_ok=True)
        _torch.save({
            "kind": "adaptive_m_selector", "M": M, "hidden": args.hidden, "log_input": log_input,
            "state_dict": model.state_dict(), "x_mean": x_mean, "x_std": x_std,
            "y_mean": y_mean, "y_std": y_std, "feature_cols": feature_cols,
        }, args.save_pt)
        print(f"wrote {args.save_pt}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sd = model.state_dict()
    with out.open("w", encoding="utf-8") as f:
        f.write(f"input_dim {M}\n")
        f.write(f"hidden {args.hidden}\n")
        f.write(f"log_input {1 if log_input else 0}\n")
        f.write(f"y_mean {repr(y_mean)}\n")
        f.write(f"y_std {repr(y_std)}\n")
        write_vec(f, "x_mean", x_mean)
        write_vec(f, "x_std", x_std)
        write_vec(f, "W1", sd["0.weight"].numpy())
        write_vec(f, "b1", sd["0.bias"].numpy())
        write_vec(f, "W2", sd["2.weight"].numpy())
        write_vec(f, "b2", sd["2.bias"].numpy())
        write_vec(f, "W3", sd["4.weight"].numpy())
        write_vec(f, "b3", sd["4.bias"].numpy())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
