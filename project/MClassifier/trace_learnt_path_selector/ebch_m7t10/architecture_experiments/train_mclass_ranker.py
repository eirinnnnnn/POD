#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np


def load_metric_dataset(path):
    path = Path(path)
    meta_path = Path(str(path) + ".meta.json")
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.shape == ():
        data = data.reshape(1)

    names = data.dtype.names
    n = int(meta.get("n", len([c for c in names if c.startswith("y_")])))
    M = int(meta.get("M", len([c for c in names if c.startswith("metric_") and c[7:].isdigit()])))
    y_cols = [f"y_{i}" for i in range(n)]
    metric_cols = [f"metric_{i}" for i in range(M)]
    X = np.stack([data[c].astype(np.float32) for c in y_cols], axis=1)
    metrics = np.stack([data[c].astype(np.float32) for c in metric_cols], axis=1)
    label = data["label"].astype(np.int64)
    return X, metrics, label, M, meta


def standardize(train_X, test_X):
    mean = train_X.mean(axis=0, keepdims=True)
    std = train_X.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return (train_X - mean) / std, (test_X - mean) / std, mean, std


def soft_targets_from_metrics(metrics, tau):
    shifted = metrics - metrics.min(axis=1, keepdims=True)
    shifted = np.clip(shifted, 0.0, 1e6)
    logits = -shifted / tau
    logits -= logits.max(axis=1, keepdims=True)
    exp_logits = np.exp(logits).astype(np.float32)
    exp_logits /= exp_logits.sum(axis=1, keepdims=True)
    return exp_logits


def main():
    parser = argparse.ArgumentParser(description="Train a metric-soft branch ranker.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--tau", type=float, default=0.25, help="Softmax temperature on metric gaps.")
    parser.add_argument("--hard-ce-weight", type=float, default=0.0, help="Optional argmin CE auxiliary weight.")
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, metrics, label, M, meta = load_metric_dataset(args.dataset)
    targets = soft_targets_from_metrics(metrics, args.tau)

    idx = np.arange(len(label))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    X_train, X_val, mean, std = standardize(X[train_idx], X[val_idx])
    t_train = targets[train_idx]
    t_val = targets[val_idx]
    y_train = label[train_idx]
    y_val = label[val_idx]

    model = nn.Sequential(
        nn.Linear(X.shape[1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, M),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    train_loader = DataLoader(
        TensorDataset(
            torch.from_numpy(X_train.astype(np.float32)),
            torch.from_numpy(t_train.astype(np.float32)),
            torch.from_numpy(y_train.astype(np.int64)),
        ),
        batch_size=args.batch,
        shuffle=True,
    )

    Xv = torch.from_numpy(X_val.astype(np.float32))
    tv = torch.from_numpy(t_val.astype(np.float32))
    yv = torch.from_numpy(y_val.astype(np.int64))
    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        hard_correct = 0
        for xb, tb, yb in train_loader:
            opt.zero_grad()
            scores = model(xb)
            log_probs = F.log_softmax(scores, dim=1)
            loss = F.kl_div(log_probs, tb, reduction="batchmean")
            if args.hard_ce_weight:
                loss = loss + args.hard_ce_weight * F.cross_entropy(scores, yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * len(yb)
            total += len(yb)
            hard_correct += int((scores.argmax(dim=1) == yb).sum().item())

        model.eval()
        with torch.no_grad():
            val_scores = model(Xv)
            val_log_probs = F.log_softmax(val_scores, dim=1)
            val_loss = float(F.kl_div(val_log_probs, tv, reduction="batchmean").item())
            val_acc = float((val_scores.argmax(dim=1) == yv).float().mean().item())
            val_top4 = float((val_scores.topk(k=min(4, M), dim=1).indices == yv[:, None]).any(dim=1).float().mean().item())
            val_top8 = float((val_scores.topk(k=min(8, M), dim=1).indices == yv[:, None]).any(dim=1).float().mean().item())

        print(
            f"epoch={epoch} train_kl={total_loss/total:.6f} "
            f"train_argmin_acc={hard_correct/total:.6f} val_kl={val_loss:.6f} "
            f"val_argmin_acc={val_acc:.6f} val_argmin_top4={val_top4:.6f} val_argmin_top8={val_top8:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "kind": "metric_soft_ranker",
            "state_dict": model.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "tau": args.tau,
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
