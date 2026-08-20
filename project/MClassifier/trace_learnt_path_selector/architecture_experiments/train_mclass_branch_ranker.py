#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_ranker import load_metric_dataset, soft_targets_from_metrics


def load_orders(path):
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=np.int64)
    if data.shape == ():
        data = data.reshape(1)
    M = len(data)
    order_cols = [f"order_{i}" for i in range(len(data.dtype.names) - 1)]
    return np.stack([data[c] for c in order_cols], axis=1).astype(np.int64)


def permute_batch(X, orders):
    # Decoder mapping is target = order[source], so permuted[target] = y[source].
    B, n = X.shape
    M = orders.shape[0]
    out = np.empty((B, M, n), dtype=np.float32)
    for b in range(M):
        out[:, b, orders[b]] = X
    return out


def main():
    parser = argparse.ArgumentParser(description="Train shared score(permuted-y) branch ranker.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--orders", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=35)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--tau", type=float, default=0.25)
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
    orders = load_orders(args.orders)
    if orders.shape[0] != M or orders.shape[1] != X.shape[1]:
        raise SystemExit(f"orders shape {orders.shape} incompatible with dataset n={X.shape[1]} M={M}")

    targets = soft_targets_from_metrics(metrics, args.tau)
    idx = np.arange(len(label))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    mean = X[train_idx].mean(axis=0, keepdims=True)
    std = X[train_idx].std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    Xn = ((X - mean) / std).astype(np.float32)

    P_train = permute_batch(Xn[train_idx], orders)
    P_val = permute_batch(Xn[val_idx], orders)
    t_train = targets[train_idx].astype(np.float32)
    t_val = targets[val_idx].astype(np.float32)
    y_train = label[train_idx].astype(np.int64)
    y_val = label[val_idx].astype(np.int64)

    scorer = nn.Sequential(
        nn.Linear(X.shape[1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    opt = torch.optim.Adam(scorer.parameters(), lr=args.lr)
    loader = DataLoader(
        TensorDataset(torch.from_numpy(P_train), torch.from_numpy(t_train), torch.from_numpy(y_train)),
        batch_size=args.batch,
        shuffle=True,
    )
    Pv = torch.from_numpy(P_val)
    tv = torch.from_numpy(t_val)
    yv = torch.from_numpy(y_val)
    for epoch in range(1, args.epochs + 1):
        scorer.train()
        total_loss = 0.0
        total = 0
        correct = 0
        for xb, tb, yb in loader:
            opt.zero_grad()
            B = xb.shape[0]
            scores = scorer(xb.reshape(B * M, -1)).reshape(B, M)
            loss = F.kl_div(F.log_softmax(scores, dim=1), tb, reduction="batchmean")
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * B
            total += B
            correct += int((scores.argmax(dim=1) == yb).sum().item())
        scorer.eval()
        with torch.no_grad():
            Bv = Pv.shape[0]
            val_scores = scorer(Pv.reshape(Bv * M, -1)).reshape(Bv, M)
            val_loss = float(F.kl_div(F.log_softmax(val_scores, dim=1), tv, reduction="batchmean").item())
            val_acc = float((val_scores.argmax(dim=1) == yv).float().mean().item())
            val_top8 = float((val_scores.topk(k=min(8, M), dim=1).indices == yv[:, None]).any(dim=1).float().mean().item())
        print(
            f"epoch={epoch} train_kl={total_loss/total:.6f} train_argmin_acc={correct/total:.6f} "
            f"val_kl={val_loss:.6f} val_argmin_acc={val_acc:.6f} val_argmin_top8={val_top8:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "kind": "shared_permuted_branch_ranker",
            "state_dict": scorer.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "orders": orders,
            "tau": args.tau,
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
