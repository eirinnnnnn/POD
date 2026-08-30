#!/usr/bin/env python3
"""model_out variant of the per-prefix-length BCE reformulation of the
adaptive-m selector: instead of regressing the scalar m_required_basin
directly (train_adaptive_m.py), predict, for every prefix length
k=1..M, the binary question "does the top-k prefix already contain a
basin hit" -- y_k = 1[k >= m_required_basin]. Derived entirely from the
existing m_required_basin column, no new data needed.

Model: same 3-layer MLP shape as train_adaptive_m.py, but the final
layer outputs M logits (one per k) instead of 1 scalar. Loss: BCE
averaged over all (sample, k) pairs. At inference, hat_m = smallest k
where sigmoid(logit_k) crosses 0.5, after a cumulative-max monotonicity
fix (nothing in the BCE loss forces the raw logits to be non-decreasing
in k).

m_required_basin is severely imbalanced (~70% of samples have m=1,
median=1, but the tail reaches M) -- confirmed to be the root cause of
both this design's earlier training instability (easy majority-k terms
dilute the boundary-relevant gradient) and the original scalar-regression
design's instability (rare large-m samples dominate squared-error
gradients when they do appear). Two mitigations added here:
  - inverse-frequency (binned) per-SAMPLE loss weighting during training,
    so all M of a sample's per-k terms get the same weight based on how
    rare its own m_required_basin value is -- keeps each per-k BCE term
    individually bounded (unlike squared error) while still fixing the
    aggregation-level dilution.
  - stratified train/val split by the same bins, so val_bce/val_miss_rate
    (used for best-epoch selection) isn't noise-dominated by which rare
    tail samples happened to land in a random 20% draw. Validation
    metrics themselves stay UNWEIGHTED (natural distribution) since
    that's what real BLER actually reflects.

Input/target come from the same CSV as train_adaptive_m.py:
    sample, M, m_required, full_ped_correct, m_required_basin,
    any_basin_exists, t{T}_trace_prob_sorted_1..M

Run with --no-log-input (default here): trace_prob is already a bounded
[0,1] sigmoid output, log-compressing it is the wrong transform, not a
neutral no-op -- same reasoning as train_adaptive_m.py's model_out
variant.
"""
import argparse
import copy
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
    m_required_basin = raw[:, col["m_required_basin"]].astype(np.float64)
    full_ped_correct = raw[:, col["any_basin_exists"]].astype(np.float32)

    k = np.arange(1, M + 1, dtype=np.float64)[None, :]
    y_prefix = (k >= m_required_basin[:, None]).astype(np.float32)  # [N, M]
    return X, y_prefix, m_required_basin, full_ped_correct, M, feature_cols


def log_transform_inputs(X):
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


def predicted_m(probs, M):
    # probs: [N, M] sigmoid outputs. Enforce monotonicity via cumulative
    # max over k, then take the smallest k crossing 0.5; if none cross,
    # clamp to M.
    cm = np.maximum.accumulate(probs, axis=1)
    crosses = cm > 0.5
    any_cross = crosses.any(axis=1)
    first_k = crosses.argmax(axis=1) + 1  # 1-indexed
    return np.where(any_cross, first_k, M).astype(np.float64)


def bin_of(y):
    # 0..9 get their own bin; y>=10 bucketed by floor(log2(y)) so the
    # long, sparse tail (singleton values up to M) still gets robust
    # per-bin frequency counts instead of each huge-m sample being its
    # own "class of one".
    y = np.asarray(y)
    return np.where(y < 10, y, 10 + np.floor(np.log2(np.maximum(y, 10)))).astype(np.int64)


def stratified_split(y, val_frac, rng):
    bins = bin_of(y)
    val_parts, train_parts = [], []
    for b in np.unique(bins):
        idx_b = np.where(bins == b)[0]
        rng.shuffle(idx_b)
        vc = 1 if len(idx_b) > 1 else 0
        vc = max(vc, int(round(len(idx_b) * val_frac)))
        vc = min(vc, len(idx_b) - 1) if len(idx_b) > 1 else 0
        val_parts.append(idx_b[:vc])
        train_parts.append(idx_b[vc:])
    val_idx = np.concatenate(val_parts) if val_parts else np.array([], dtype=np.int64)
    train_idx = np.concatenate(train_parts) if train_parts else np.array([], dtype=np.int64)
    rng.shuffle(val_idx)
    rng.shuffle(train_idx)
    return train_idx, val_idx


def sample_weights(y_train):
    bins = bin_of(y_train)
    vals, counts = np.unique(bins, return_counts=True)
    freq = dict(zip(vals.tolist(), counts.tolist()))
    w = np.array([1.0 / freq[b] for b in bins], dtype=np.float64)
    w = w / w.mean()  # normalize so average weight is 1 (keeps loss scale comparable)
    return w.astype(np.float32)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True)
    parser.add_argument("--out", required=True, help="plain-text weight file path")
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--log-input", action="store_true",
                         help="enable the log1p input transform (off by default for model_out -- "
                              "trace_prob is already bounded [0,1])")
    parser.add_argument("--no-weighting", action="store_true",
                         help="disable inverse-frequency sample weighting and stratified split "
                              "(for A/B comparison against the balanced version)")
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, y_prefix, m_req_basin, full_ped_correct, M, feature_cols = load_dataset(args.data)
    log_input = args.log_input
    if log_input:
        X = log_transform_inputs(X)
    n = X.shape[0]

    if args.no_weighting:
        idx = np.arange(n)
        rng.shuffle(idx)
        val_count = max(1, int(n * args.val_frac))
        val_idx, train_idx = idx[:val_count], idx[val_count:]
    else:
        train_idx, val_idx = stratified_split(m_req_basin, args.val_frac, rng)

    X_train, X_val, x_mean, x_std = standardize(X[train_idx], X[val_idx])
    y_train = y_prefix[train_idx]
    y_val = y_prefix[val_idx]
    m_req_val = m_req_basin[val_idx]
    fpc_val = full_ped_correct[val_idx]

    if args.no_weighting:
        weights_train = np.ones(len(train_idx), dtype=np.float32)
    else:
        weights_train = sample_weights(m_req_basin[train_idx])

    model = nn.Sequential(
        nn.Linear(M, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, M),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    loader = DataLoader(
        TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train), torch.from_numpy(weights_train)),
        batch_size=args.batch_size, shuffle=True,
    )
    Xv = torch.from_numpy(X_val)
    yv = torch.from_numpy(y_val)
    bce_elem = torch.nn.BCEWithLogitsLoss(reduction="none")
    bce = torch.nn.BCEWithLogitsLoss()

    best_val_loss = float("inf")
    best_epoch = 0
    best_state = None

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss, total = 0.0, 0
        for xb, yb, wb in loader:
            opt.zero_grad()
            logits = model(xb)
            per_sample = bce_elem(logits, yb).mean(dim=1)  # [B]
            loss = (per_sample * wb).sum() / wb.sum()
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * xb.shape[0]
            total += xb.shape[0]

        model.eval()
        with torch.no_grad():
            logits_v = model(Xv)
            val_loss = float(bce(logits_v, yv).item())  # unweighted -- natural distribution
            probs_v = torch.sigmoid(logits_v).numpy()
            used_m = predicted_m(probs_v, M)
            miss = float((((used_m < m_req_val) & (fpc_val > 0.5))).mean())
            mean_m_used = float(used_m.mean())
        print(f"epoch={epoch} train_loss={total_loss/total:.6f} val_bce={val_loss:.6f} "
              f"val_miss_rate={miss:.4f} val_mean_m={mean_m_used:.2f} (M={M})")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())

    print(f"selecting best epoch={best_epoch} (val_bce={best_val_loss:.6f}) for export")
    model.load_state_dict(best_state)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sd = model.state_dict()
    with out.open("w", encoding="utf-8") as f:
        f.write(f"input_dim {M}\n")
        f.write(f"hidden {args.hidden}\n")
        f.write(f"output_dim {M}\n")
        f.write(f"log_input {1 if log_input else 0}\n")
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
