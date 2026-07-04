#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_branch_ranker import load_orders
from train_mclass_ranker import load_metric_dataset


def target_scale_from_gaps(gaps):
    positive = gaps[gaps > 1e-9]
    if positive.size == 0:
        return 1.0
    scale = float(np.median(positive))
    return scale if scale > 1e-9 else 1.0


def make_pair_features(X, orders, branch_idx):
    order = orders[branch_idx]
    permuted = np.empty_like(X)
    permuted[:, order] = X
    order_norm = order.astype(np.float32) / max(1, X.shape[1] - 1)
    order_features = np.broadcast_to(order_norm.reshape(1, -1), X.shape)
    return np.concatenate([X, permuted, order_features], axis=1).astype(np.float32)


def main():
    parser = argparse.ArgumentParser(description="Train g(y, permutation) -> final metric gap.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--orders", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--hidden", type=int, default=256)
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=512)
    parser.add_argument("--val-frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--branches-per-sample", type=int, default=16)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch.utils.data import DataLoader, Dataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, metrics, label, M, meta = load_metric_dataset(args.dataset)
    orders = load_orders(args.orders)
    if orders.shape[0] != M or orders.shape[1] != X.shape[1]:
        raise SystemExit(f"orders shape {orders.shape} incompatible with dataset n={X.shape[1]} M={M}")

    idx = np.arange(len(label))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    mean = X[train_idx].mean(axis=0, keepdims=True)
    std = X[train_idx].std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    Xn = ((X - mean) / std).astype(np.float32)

    gaps = metrics - metrics.min(axis=1, keepdims=True)
    scale = target_scale_from_gaps(gaps[train_idx])
    targets = np.log1p(gaps / scale).astype(np.float32)

    class PairDataset(Dataset):
        def __init__(self, sample_idx, random_branches):
            self.sample_idx = np.asarray(sample_idx, dtype=np.int64)
            self.random_branches = random_branches
            if random_branches:
                self.length = len(self.sample_idx) * args.branches_per_sample
            else:
                self.length = len(self.sample_idx) * M

        def __len__(self):
            return self.length

        def __getitem__(self, i):
            s = self.sample_idx[i // (args.branches_per_sample if self.random_branches else M)]
            if self.random_branches:
                # Always include the true argmin in one slot on average by cycling slot 0.
                slot = i % args.branches_per_sample
                b = int(label[s]) if slot == 0 else int(rng.integers(0, M))
            else:
                b = i % M
            y = Xn[s]
            py = np.empty_like(y)
            py[orders[b]] = y
            order_norm = orders[b].astype(np.float32) / max(1, X.shape[1] - 1)
            feat = np.concatenate([y, py, order_norm]).astype(np.float32)
            return feat, np.float32(targets[s, b]), np.int64(s), np.int64(b)

    model = nn.Sequential(
        nn.Linear(X.shape[1] * 3, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    train_loader = DataLoader(
        PairDataset(train_idx, random_branches=True),
        batch_size=args.batch,
        shuffle=True,
    )

    def evaluate_argmin(sample_idx):
        model.eval()
        rows = []
        with torch.no_grad():
            for b in range(M):
                feat = make_pair_features(Xn[sample_idx], orders, b)
                pred = model(torch.from_numpy(feat)).squeeze(1).cpu().numpy()
                rows.append(pred)
        pred_scores = np.stack(rows, axis=1)
        pred_argmin = pred_scores.argmin(axis=1)
        true_argmin = label[sample_idx]
        top8 = np.argsort(pred_scores, axis=1)[:, :min(8, M)]
        return (
            float((pred_argmin == true_argmin).mean()),
            float((top8 == true_argmin[:, None]).any(axis=1).mean()),
        )

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        for xb, tb, _, _ in train_loader:
            opt.zero_grad()
            pred = model(xb).squeeze(1)
            loss = F.smooth_l1_loss(pred, tb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * len(tb)
            total += len(tb)
        val_acc, val_top8 = evaluate_argmin(val_idx)
        print(
            f"epoch={epoch} train_huber={total_loss/total:.6f} "
            f"val_argmin_acc={val_acc:.6f} val_argmin_top8={val_top8:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    import torch
    torch.save(
        {
            "kind": "explicit_y_permutation_metric_gap_regressor",
            "state_dict": model.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "orders": orders,
            "target_scale": scale,
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
