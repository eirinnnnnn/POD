#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_pair_ranker import make_pair_features
from train_mclass_ranker import load_metric_dataset


def main():
    parser = argparse.ArgumentParser(description="Predict top-k branches with explicit (y, permutation) regressor.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--topk", type=int, default=1)
    parser.add_argument("--batch", type=int, default=256)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
    except ImportError as exc:
        raise SystemExit("PyTorch is required for prediction: pip install torch") from exc

    X, _, _, _, _ = load_metric_dataset(args.dataset)
    ckpt = torch.load(args.model, map_location="cpu")
    n = int(ckpt["n"])
    M = int(ckpt["M"])
    orders = np.asarray(ckpt["orders"], dtype=np.int64)
    model = nn.Sequential(
        nn.Linear(n * 3, int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), 1),
    )
    model.load_state_dict(ckpt["state_dict"])
    model.eval()

    Xn = ((X - ckpt["mean"]) / ckpt["std"]).astype(np.float32)
    scores = np.empty((len(Xn), M), dtype=np.float32)
    with torch.no_grad():
        for b in range(M):
            for start in range(0, len(Xn), args.batch):
                xb = Xn[start:start + args.batch]
                feat = make_pair_features(xb, orders, b)
                pred = model(torch.from_numpy(feat)).squeeze(1).cpu().numpy()
                scores[start:start + len(xb), b] = pred
    topk = np.argsort(scores, axis=1)[:, :min(args.topk, M)]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("sample")
        for i in range(topk.shape[1]):
            f.write(f",idx{i}")
        f.write("\n")
        for s, row in enumerate(topk):
            f.write(str(s))
            for idx in row:
                f.write(f",{int(idx)}")
            f.write("\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
