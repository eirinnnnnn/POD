#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_ranker import load_metric_dataset
from train_mclass_branch_ranker import permute_batch


def main():
    parser = argparse.ArgumentParser(description="Predict top-k branches with shared permuted-y scorer.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--topk", type=int, default=1)
    parser.add_argument("--batch", type=int, default=128)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
    except ImportError as exc:
        raise SystemExit("PyTorch is required for prediction: pip install torch") from exc

    X, _, _, _, _ = load_metric_dataset(args.dataset)
    ckpt = torch.load(args.model, map_location="cpu")
    M = int(ckpt["M"])
    n = int(ckpt["n"])
    orders = np.asarray(ckpt["orders"], dtype=np.int64)
    scorer = nn.Sequential(
        nn.Linear(n, int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), 1),
    )
    scorer.load_state_dict(ckpt["state_dict"])
    scorer.eval()

    Xn = ((X - ckpt["mean"]) / ckpt["std"]).astype(np.float32)
    top_rows = []
    with torch.no_grad():
        for start in range(0, len(Xn), args.batch):
            xb = Xn[start:start + args.batch]
            pb = permute_batch(xb, orders)
            B = pb.shape[0]
            scores = scorer(torch.from_numpy(pb).reshape(B * M, n)).reshape(B, M)
            top = torch.topk(scores, k=min(args.topk, M), dim=1).indices.cpu().numpy()
            top_rows.append(top)
    topk = np.concatenate(top_rows, axis=0)

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
