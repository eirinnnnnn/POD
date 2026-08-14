#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass import standardize
from train_mclass_branch_ranker import load_orders, permute_batch
from train_mclass_ranker import load_metric_dataset


def synthetic_labels(X, orders, window):
    # Decoder mapping is target = order[source], so permuted[target] = y[source].
    P = permute_batch(X.astype(np.float32), orders)
    w = min(window, X.shape[1])
    scores = np.mean(np.abs(P[:, :, :w]), axis=2)
    return scores.argmax(axis=1).astype(np.int64), scores.astype(np.float32)


def topk_acc(scores, y, ks):
    out = {}
    for k in ks:
        kk = min(k, scores.shape[1])
        top = np.argpartition(-scores, kk - 1, axis=1)[:, :kk]
        out[k] = float((top == y[:, None]).any(axis=1).mean())
    return out


def main():
    parser = argparse.ArgumentParser(description="Synthetic branch-rule learnability check.")
    parser.add_argument("--train", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--orders", required=True)
    parser.add_argument("--mode", choices=["plain", "branch"], default="plain")
    parser.add_argument("--branch-loss", choices=["ce", "mse"], default="ce")
    parser.add_argument("--out", required=True)
    parser.add_argument("--window", type=int, default=16)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X_train_raw, _, _, M, _ = load_metric_dataset(args.train)
    X_test_raw, _, _, M_test, _ = load_metric_dataset(args.test)
    orders = load_orders(args.orders)
    if M != M_test or orders.shape != (M, X_train_raw.shape[1]):
        raise SystemExit("dataset/order shape mismatch")

    y_train, target_scores_train = synthetic_labels(X_train_raw, orders, args.window)
    y_test, target_scores_test = synthetic_labels(X_test_raw, orders, args.window)

    X_train, X_test, mean, std = standardize(X_train_raw, X_test_raw)
    X_train = X_train.astype(np.float32)
    X_test = X_test.astype(np.float32)

    idx = np.arange(len(y_train))
    rng.shuffle(idx)
    split = max(1, int(0.8 * len(idx)))
    train_idx = idx[:split]
    val_idx = idx[split:]

    if args.mode == "plain":
        model = nn.Sequential(
            nn.Linear(X_train.shape[1], args.hidden),
            nn.ReLU(),
            nn.Linear(args.hidden, args.hidden),
            nn.ReLU(),
            nn.Linear(args.hidden, M),
        )
        loader = DataLoader(
            TensorDataset(torch.from_numpy(X_train[train_idx]), torch.from_numpy(y_train[train_idx])),
            batch_size=args.batch,
            shuffle=True,
        )

        def score_numpy(X):
            model.eval()
            rows = []
            with torch.no_grad():
                for start in range(0, len(X), 1024):
                    rows.append(model(torch.from_numpy(X[start:start + 1024])).cpu().numpy())
            return np.concatenate(rows, axis=0)

        opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
        loss_name = "ce"
        for epoch in range(1, args.epochs + 1):
            model.train()
            total_loss = 0.0
            total = 0
            for xb, yb in loader:
                opt.zero_grad()
                logits = model(xb)
                loss = F.cross_entropy(logits, yb)
                loss.backward()
                opt.step()
                total_loss += float(loss.item()) * len(yb)
                total += len(yb)

            val_scores = score_numpy(X_train[val_idx])
            test_scores = score_numpy(X_test)
            val_top = topk_acc(val_scores, y_train[val_idx], [1, 4, 8])
            test_top = topk_acc(test_scores, y_test, [1, 4, 8])
            print(
                f"epoch={epoch} train_{loss_name}={total_loss/total:.6f} "
                f"val_top1={val_top[1]:.6f} val_top4={val_top[4]:.6f} val_top8={val_top[8]:.6f} "
                f"test_top1={test_top[1]:.6f} test_top4={test_top[4]:.6f} test_top8={test_top[8]:.6f}"
            )
    else:
        P_train = permute_batch(X_train, orders)
        P_test = permute_batch(X_test, orders)
        model = nn.Sequential(
            nn.Linear(X_train.shape[1], args.hidden),
            nn.ReLU(),
            nn.Linear(args.hidden, args.hidden),
            nn.ReLU(),
            nn.Linear(args.hidden, 1),
        )
        loader = DataLoader(
            TensorDataset(
                torch.from_numpy(P_train[train_idx]),
                torch.from_numpy(y_train[train_idx]),
                torch.from_numpy(target_scores_train[train_idx]),
            ),
            batch_size=max(16, args.batch // 4),
            shuffle=True,
        )

        def score_permuted_numpy(P):
            model.eval()
            rows = []
            with torch.no_grad():
                for start in range(0, len(P), 256):
                    xb = torch.from_numpy(P[start:start + 256])
                    B = xb.shape[0]
                    rows.append(model(xb.reshape(B * M, -1)).reshape(B, M).cpu().numpy())
            return np.concatenate(rows, axis=0)

        opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
        for epoch in range(1, args.epochs + 1):
            model.train()
            total_loss = 0.0
            total = 0
            for xb, yb, sb in loader:
                opt.zero_grad()
                B = xb.shape[0]
                logits = model(xb.reshape(B * M, -1)).reshape(B, M)
                if args.branch_loss == "mse":
                    loss = F.mse_loss(logits, sb)
                else:
                    loss = F.cross_entropy(logits, yb)
                loss.backward()
                opt.step()
                total_loss += float(loss.item()) * B
                total += B

            val_scores = score_permuted_numpy(P_train[val_idx])
            test_scores = score_permuted_numpy(P_test)
            val_top = topk_acc(val_scores, y_train[val_idx], [1, 4, 8])
            test_top = topk_acc(test_scores, y_test, [1, 4, 8])
            print(
                f"epoch={epoch} train_{args.branch_loss}={total_loss/total:.6f} "
                f"val_top1={val_top[1]:.6f} val_top4={val_top[4]:.6f} val_top8={val_top[8]:.6f} "
                f"test_top1={test_top[1]:.6f} test_top4={test_top[4]:.6f} test_top8={test_top[8]:.6f}"
            )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "kind": "synthetic_branch_rule",
            "mode": args.mode,
            "window": args.window,
            "state_dict": model.state_dict(),
            "n": X_train.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "orders": orders,
        },
        out,
    )
    print(f"wrote {out}")
    print("synthetic_label_top10_train", np.bincount(y_train, minlength=M).argsort()[-10:][::-1].tolist())
    print("synthetic_label_top10_test", np.bincount(y_test, minlength=M).argsort()[-10:][::-1].tolist())


if __name__ == "__main__":
    main()
