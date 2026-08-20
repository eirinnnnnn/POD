#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np


def load_dataset(path):
    path = Path(path)
    meta_path = Path(str(path) + ".meta.json")
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.shape == ():
        data = data.reshape(1)

    n = int(meta.get("n", len([c for c in data.dtype.names if c.startswith("y_")])))
    y_cols = [f"y_{i}" for i in range(n)]
    X = np.stack([data[c].astype(np.float32) for c in y_cols], axis=1)
    y = data["label"].astype(np.int64)
    M = int(meta.get("M", int(y.max()) + 1))
    return X, y, M, meta


def standardize(train_X, test_X):
    mean = train_X.mean(axis=0, keepdims=True)
    std = train_X.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return (train_X - mean) / std, (test_X - mean) / std, mean, std


def main():
    parser = argparse.ArgumentParser(description="Train an M-classifier MLP.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", required=True, help="Output .npz model path.")
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, y, M, meta = load_dataset(args.dataset)
    idx = np.arange(len(y))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    X_train, X_val, mean, std = standardize(X[train_idx], X[val_idx])
    y_train = y[train_idx]
    y_val = y[val_idx]

    model = nn.Sequential(
        nn.Linear(X.shape[1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, M),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    loss_fn = nn.CrossEntropyLoss()

    train_loader = DataLoader(
        TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train)),
        batch_size=args.batch,
        shuffle=True,
    )

    Xv = torch.from_numpy(X_val)
    yv = torch.from_numpy(y_val)
    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        correct = 0
        for xb, yb in train_loader:
            opt.zero_grad()
            logits = model(xb)
            loss = loss_fn(logits, yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * len(yb)
            total += len(yb)
            correct += int((logits.argmax(dim=1) == yb).sum().item())

        model.eval()
        with torch.no_grad():
            val_logits = model(Xv)
            val_loss = float(loss_fn(val_logits, yv).item())
            val_acc = float((val_logits.argmax(dim=1) == yv).float().mean().item())

        print(
            f"epoch={epoch} train_loss={total_loss/total:.6f} "
            f"train_acc={correct/total:.6f} val_loss={val_loss:.6f} val_acc={val_acc:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dict": model.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
