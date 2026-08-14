#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_trace_ranker import load_trace_dataset


def main():
    parser = argparse.ArgumentParser(description="Predict top-k branches from decoder trace features.")
    parser.add_argument("--trace", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--topk", type=int, default=64)
    parser.add_argument("--batch-samples", type=int, default=256)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
    except ImportError as exc:
        raise SystemExit("PyTorch is required for prediction: pip install torch") from exc

    ckpt = torch.load(args.model, map_location="cpu", weights_only=False)
    X, _, _, _, M, _ = load_trace_dataset(args.trace, branch_onehot=bool(ckpt["branch_onehot"]))
    if M != int(ckpt["M"]):
        raise SystemExit(f"trace M={M} but model M={ckpt['M']}")

    model = nn.Sequential(
        nn.Linear(int(ckpt["input_dim"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), 1),
    )
    model.load_state_dict(ckpt["state_dict"])
    model.eval()

    mean = np.asarray(ckpt["mean"], dtype=np.float32).reshape(1, 1, -1)
    std = np.asarray(ckpt["std"], dtype=np.float32).reshape(1, 1, -1)
    Xn = ((X - mean) / std).astype(np.float32)
    largest = ckpt["direction"] == "largest"

    top_rows = []
    with torch.no_grad():
        for start in range(0, Xn.shape[0], args.batch_samples):
            xb = torch.from_numpy(Xn[start:start + args.batch_samples])
            B = xb.shape[0]
            scores = model(xb.reshape(B * M, -1)).reshape(B, M)
            top = torch.topk(scores, k=min(args.topk, M), dim=1, largest=largest).indices.cpu().numpy()
            top_rows.append(top)
    topk = np.concatenate(top_rows, axis=0)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("sample")
        for i in range(topk.shape[1]):
            f.write(f",idx{i}")
        f.write("\n")
        for sample, row in enumerate(topk):
            f.write(str(sample))
            for idx in row:
                f.write(f",{int(idx)}")
            f.write("\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
