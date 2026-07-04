#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np


TARGET_COLS = {"sample", "branch", "metric_gap", "basin", "correct", "final_metric"}


def load_trace_dataset(path, branch_onehot=True):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=np.float32)
    names = list(data.dtype.names)
    samples = data["sample"].astype(np.int64)
    branches = data["branch"].astype(np.int64)
    n_samples = int(samples.max()) + 1
    M = int(branches.max()) + 1
    if len(samples) != n_samples * M:
        raise ValueError(f"expected dense sample/branch rows, got {len(samples)} for {n_samples}*{M}")

    order = np.lexsort((branches, samples))
    data = data[order]
    branches = data["branch"].astype(np.int64).reshape(n_samples, M)
    if not np.all(branches == np.arange(M)[None, :]):
        raise ValueError("trace rows are not dense in branch order")

    feature_cols = [name for name in names if name not in TARGET_COLS]
    X = np.stack([data[name] for name in feature_cols], axis=1).reshape(n_samples, M, -1)
    if branch_onehot:
        eye = np.eye(M, dtype=np.float32)[np.arange(M)]
        eye = np.broadcast_to(eye[None, :, :], (n_samples, M, M))
        X = np.concatenate([X.astype(np.float32), eye], axis=2)
        feature_cols = feature_cols + [f"branch_{i}" for i in range(M)]
    gaps = data["metric_gap"].reshape(n_samples, M).astype(np.float32)
    basin = data["basin"].reshape(n_samples, M).astype(np.float32)
    labels = gaps.argmin(axis=1).astype(np.int64)
    return X.astype(np.float32), gaps, basin, labels, M, feature_cols


def standardize_rows(train, val):
    flat = train.reshape(-1, train.shape[-1])
    mean = flat.mean(axis=0, keepdims=True)
    std = flat.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return (train - mean) / std, (val - mean) / std, mean.reshape(-1), std.reshape(-1)


def gap_scale(gaps):
    positive = gaps[gaps > 1e-9]
    if positive.size == 0:
        return 1.0
    return max(float(np.median(positive)), 1e-9)


def main():
    parser = argparse.ArgumentParser(description="Train a branch scorer from decoder trace features.")
    parser.add_argument("--trace", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--target", choices=["basin", "log_gap"], default="basin")
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch-rows", type=int, default=4096)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--no-branch-onehot", action="store_true")
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, gaps, basin, labels, M, feature_cols = load_trace_dataset(
        args.trace, branch_onehot=not args.no_branch_onehot
    )
    idx = np.arange(X.shape[0])
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    X_train, X_val, mean, std = standardize_rows(X[train_idx], X[val_idx])
    if args.target == "basin":
        y_train = basin[train_idx]
        y_val = basin[val_idx]
        direction = "largest"
        scale = 1.0
    else:
        scale = gap_scale(gaps[train_idx])
        y_train = np.log1p(gaps[train_idx] / scale).astype(np.float32)
        y_val = np.log1p(gaps[val_idx] / scale).astype(np.float32)
        direction = "smallest"

    model = nn.Sequential(
        nn.Linear(X.shape[-1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    loader = DataLoader(
        TensorDataset(
            torch.from_numpy(X_train.reshape(-1, X.shape[-1]).astype(np.float32)),
            torch.from_numpy(y_train.reshape(-1).astype(np.float32)),
        ),
        batch_size=args.batch_rows,
        shuffle=True,
    )

    Xv = torch.from_numpy(X_val.astype(np.float32))
    yv = torch.from_numpy(y_val.astype(np.float32))
    label_v = torch.from_numpy(labels[val_idx].astype(np.int64))
    basin_v = torch.from_numpy(basin[val_idx].astype(np.float32))

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        for xb, yb in loader:
            opt.zero_grad()
            pred = model(xb).reshape(-1)
            if args.target == "basin":
                loss = F.binary_cross_entropy_with_logits(pred, yb)
            else:
                loss = F.smooth_l1_loss(pred, yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * xb.shape[0]
            total += xb.shape[0]

        model.eval()
        with torch.no_grad():
            Bv = Xv.shape[0]
            scores = model(Xv.reshape(Bv * M, -1)).reshape(Bv, M)
            if args.target == "basin":
                val_loss = float(F.binary_cross_entropy_with_logits(scores, yv).item())
                ranked = torch.topk(scores, k=M, dim=1, largest=True).indices
            else:
                val_loss = float(F.smooth_l1_loss(scores, yv).item())
                ranked = torch.topk(scores, k=M, dim=1, largest=False).indices
            acc1 = float((ranked[:, 0] == label_v).float().mean().item())
            hit1 = float((basin_v.gather(1, ranked[:, :1]).amax(dim=1) > 0.5).float().mean().item())
            hit4 = float((basin_v.gather(1, ranked[:, : min(4, M)]).amax(dim=1) > 0.5).float().mean().item())
            hit8 = float((basin_v.gather(1, ranked[:, : min(8, M)]).amax(dim=1) > 0.5).float().mean().item())

        print(
            f"epoch={epoch} train_loss={total_loss/total:.6f} val_loss={val_loss:.6f} "
            f"val_label_top1={acc1:.6f} val_basin_hit1={hit1:.6f} "
            f"val_basin_hit4={hit4:.6f} val_basin_hit8={hit8:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    import torch

    torch.save(
        {
            "kind": "trace_branch_ranker",
            "target": args.target,
            "direction": direction,
            "state_dict": model.state_dict(),
            "M": M,
            "input_dim": X.shape[-1],
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "branch_onehot": not args.no_branch_onehot,
            "feature_cols": feature_cols,
            "target_scale": scale,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
