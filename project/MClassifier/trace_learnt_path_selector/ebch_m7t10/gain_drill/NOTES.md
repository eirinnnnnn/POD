# Gain drill: make trace_learned actually beat static_pm_rank at t=64

## Goal

At checkpoint t=64 (of n=128), k=8, the shipped trace-learned model
(`results/clean_trace_upto_t64_basin_h64_e20.pt`) does NOT show a clear
win over the model-free static-pm-rank baseline -- across the whole SNR
sweep it's a coin flip, sometimes marginally better, sometimes worse,
well inside noise. This drill's job: find a genuine, reproducible gain at
t=64 specifically, so the model's value-add can actually be claimed
there (not just at earlier checkpoints where the gap is already clean).

## Baseline to beat (results/ebn0_sweep/results.txt, t=64, k=8)

| SNR  | static_pm_rank BLER | trace_learned BLER | trace wins? |
|------|---------------------|---------------------|-------------|
| 2.00 | 0.137477            | 0.139336            | no          |
| 2.25 | 0.0910091           | 0.0906424           | marginal    |
| 2.50 | 0.0534599           | 0.0542806           | no          |
| 2.75 | 0.0277347           | 0.0278912           | no          |
| 3.00 | 0.01386             | 0.01409             | no          |
| 3.50 | 0.00303             | 0.00287             | marginal    |
| 4.00 | 0.00030648          | 0.000326468         | no          |
| 4.50 | 0.0000249667        | 0.0000274634        | no          |

Trace-learned "wins" on 2/8 SNR points, and only marginally. Not a
claimable gain as-is.

## Why this is plausibly hard at t=64 specifically

`results/clean_trace_12k_upto_t64.csv` already has CUMULATIVE features
(t4, t8, t16, t32, t64 all present as columns) -- the shipped model
already sees the full mid-decode history, not just a single checkpoint,
and still doesn't beat raw pm_min. That points away from "needs more
raw signal" and toward one of:

1. **Metric saturation**: by t=64 (halfway through n=128), many AED
   branches that converge to the correct codeword have already locked
   onto bit-identical metrics (the same phenomenon documented earlier
   this session for Golay24's checkpoint sweep) -- the raw pm_min
   ordering may already be close to information-theoretically optimal
   by this point, leaving little room for a learned re-ranking to add
   value on TOP of it.
2. **Wrong training objective**: the shipped model was trained with
   `--target basin` (binary per-branch correctness classification via
   BCE). That optimizes classification accuracy, not top-k RANKING
   quality, which is what static_bler(k) actually measures.
   `architectures/trace/train_mclass_trace_ranker.py` already supports
   `--target log_gap` (regress log1p(metric_gap), a genuinely
   rank-oriented target) -- never tried at t=64. Cheapest lever to test
   first: zero new code, just a different training run.
3. **Single-SNR training mismatch**: trained only on SNR=2.0dB data
   (`data/m7t10_m64_snr2p0_train_metric_12k.csv`-equivalent), evaluated
   across 2.00-4.50dB -- the exact distribution-shift caveat flagged in
   EXPERIMENT_NOTES.txt, never resolved for the trace models.
4. **Undertrained**: only 20 epochs, hidden=64 -- both arbitrary
   defaults never swept for this specific checkpoint.

## Experiment log

**#1 -- log_gap target vs basin target** (same hidden=64/epochs=20/data):
`gain_drill/t64_log_gap_h64_e20.pt`, final epoch val_basin_hit8=0.8813
vs the shipped basin-target model's val_basin_hit8=0.8829 (both from
their own training logs, same held-out split size). Essentially
identical -- **rules out "wrong training objective" as the bottleneck**,
doesn't confirm it. Both losses converge to the same ranking quality
ceiling on this data. Next: rule out undertrained/too-small-capacity
(more epochs, wider hidden) before concluding it's a genuine information
ceiling at t=64.

**#2 -- capacity/undertraining check**: `gain_drill/t64_basin_h128_e40.pt`
(hidden=128, epochs=40, basin target) -- final val_basin_hit8=0.8808,
same ceiling as h64/e20 (0.8829) and log_gap (0.8813). Both loss curves
plateau almost immediately (epoch 1 already ~0.88) rather than climbing
over training -- the signature of a saturated model, not an undertrained
one. **Rules out capacity/epochs as the bottleneck too.**

**Frozen-bit structural check** (added `--dump-frozen` to
`adaptive_path_selector/ebch_m7t10/mclass_adaptive_m_dataset.cpp`, reusing
`MClassDecoder::relationSizes()`/`informationPositions()` which already
existed in the shared `common/mclass_decoder_wrapper.h` -- no decoder
changes needed; dump saved as `frozen_bit_dump_snr3p0.txt`):
checkpoint t is captured right after decode_idx=t-1, i.e. right BEFORE
decode_idx=t itself. All four trace-ranker checkpoints turn out to share
the exact same structural property -- each sits right before a
frozen/dynamic-frozen bit: t=8→decode_idx=8 (relation_span=2),
t=16→16 (span=4), t=32→32 (span=5), t=64→64 (span=9). So "right before a
consistency check" does NOT single out t=64 -- it's true of every
checkpoint that DOES show a real gap too. The one real structural
difference: the constraint right after t=64 is markedly WIDER
(relation_span=9 vs max 5 at the earlier checkpoints), part of an entire
frozen cluster at decode_idx 64-70 (spans 9,1,1,2,2,3,6) -- plausibly one
dominant consistency check overwhelms whatever finer-grained ranking
signal a learned model could add at t=64 specifically. Also checked:
t=64 is NOT an information-saturation point the way Golay24's late
checkpoints were -- 35 of 64 info bits are still unresolved beyond
decode_idx=64, so remaining uncertainty isn't the issue either.

**Net so far**: two architecture/objective levers ruled out; the
checkpoint's own local structure (wide constraint immediately following)
looks like the more likely lead, but isn't yet a confirmed fix. Untested:
multi-SNR training data (the one remaining hypothesis from the original
list), and deliberately probing checkpoints near t=64 (e.g. t=63, t=71 --
right before vs. right after the whole 64-70 cluster resolves) the same
way the Golay t=20/21 gap was isolated.

**#3 -- static_pm_rank-only probe around t=60-72** (SNR=3.0, k=8, no
model, mirroring the Golay t=20/21 diagnostic): static_bler(t) at
t=60,62,63,64,65,66,67,68,69,70,71,72 = 0.0144, 0.0144, 0.0144, 0.0135,
0.0138, 0.0135, 0.0138, 0.0118, 0.0123, 0.0113, 0.0095, 0.0095 (4-8k
samples per point, noisy). **No sharp cliff like Golay's t=20/21** --
instead a gradual decline, roughly flat through t~60-67 then a smoother
(not binary) drop from t~68 through t~71-72, loosely tracking where the
decode_idx=64-70 frozen cluster finishes resolving. The "checkpoint sits
right before one dominant consistency check" story that explained
Golay24's cliff does not transfer cleanly to eBCH -- the effect here is
diffuse across several decode positions, not concentrated in one jump.
This means picking a different exact checkpoint near 64 (e.g. t=71
instead of t=64) probably won't produce a Golay-style dramatic
before/after split for the trace-learned gap specifically -- it would
mostly just shift both static and trace to a slightly lower absolute
BLER together, likely still tied.

**#4 -- bracketing where the gap disappears, t=32 (known gap) to t=64
(known tie)**: trained fresh models at t=48 and t=40 (cumulative
checkpoints 4,8,16,32,{48 or 40}), same basin/hidden=64/epochs=20 recipe,
using the existing `results/clean_train_metric_12k.csv`/
`clean_test_metric_3k.csv` base datasets re-extracted at the new
checkpoint via `mclass_trace_dataset` (no resimulation needed). BLER
comparison (static vs trace, k=8, SNR=2.0/2.5/3.0, 6000 samples,
channel-seed 777/message-seed 888):

| t  | SNR=2.0 static/trace | SNR=2.5 static/trace | SNR=3.0 static/trace | gap? |
|----|----------------------|-----------------------|-----------------------|------|
| 32 | 0.3223 / 0.3049 (shipped) | --                | --                    | yes, ~5% |
| 40 | 0.2637 / 0.2592       | 0.1328 / 0.1258       | 0.0493 / 0.0460       | yes, ~2-7% |
| 48 | 0.2008 / 0.2005       | 0.0932 / 0.0900       | 0.0273 / 0.0277       | no (tied) |
| 64 | 0.1375 / 0.1393       | 0.0535 / 0.0543       | 0.0139 / 0.0141       | no (tied) |

**The real, consistent, model-favoring gap survives to t=40 and is gone
by t=48** -- narrower window than the original t=32-to-64 sweep implied.
Models/data for both live in this directory
(`clean_trace_{12k,3k}_upto_t{40,48}.csv`, `t{40,48}_basin_h64_e20.pt`,
`t{40,48}_weights.txt`). Not yet narrowed further (e.g. t=44) -- next
step if pinning down the exact transition point further is wanted.

**#5 -- PED27-SCL4 complexity reference for the t=40 result**: created
`results/ebn0_sweep_m27.ini` (aed_L=27, list_size=4, same convention as
the existing ebn0_sweep_m8/m24/m32.ini family) and decoded all 27
branches in full (no checkpoint pruning) at the same SNR points/seeds as
the t=40 comparison. Rationale: decoding all M=64 branches through
checkpoint t=40 then pruning to k=8 for the remaining 88 positions costs
about 64x40 + 8x88 = 3264 position-units, roughly matching a full decode
of M=27 branches (27x128 = 3456) -- a genuinely smaller, unpruned
ensemble at comparable total complexity.

| SNR | static(t=40,k=8) | trace(t=40,k=8) | PED27-SCL4 (full, no pruning) |
|-----|-------------------|-------------------|--------------------------------|
| 2.0 | 0.2637            | 0.2592            | **0.1845**                     |
| 2.5 | 0.1328            | 0.1258            | **0.0840**                     |
| 3.0 | 0.0493            | 0.0460            | **0.0277**                     |

PED27 beats BOTH pruning methods substantially (~29-40% lower BLER than
trace_learned). So while trace_learned's win over static_pm_rank at
t=40 is real (confirmed in #4), both remain clearly dominated by a
complexity-matched smaller-but-unpruned ensemble. The "gain" claim holds
specifically against static_pm_rank, not against this alternative -- worth
being explicit about which baseline any claim is measured against.

**#6 -- does an earlier, cheaper checkpoint change the PED verdict?**
Complexity-matched PED_M shrinks a lot with smaller t (overhead of
running all 64 branches through the checkpoint dominates the cost):
PED_M(t,k=8) = (64t + 8(128-t))/128 -- t=40 -> M~25.5, t=8 -> M~11.5.
The shipped notes' existing `full_ped_m8` comparison at t=8 (M=8, NOT
complexity-matched -- weaker/easier reference than the true M~11.5)
already showed trace < static < full_ped_m8. Built the properly matched
`results/ebn0_sweep_m12.ini` (aed_L=12) and reran at t=8/k=8 using the
already-shipped `clean_trace_upto_t8_weights.txt` (no retraining needed):

| SNR | static(t=8,k=8) | trace(t=8,k=8) | PED12-SCL4 (full, properly matched) |
|-----|-------------------|-------------------|--------------------------------------|
| 2.0 | 0.4500            | 0.3928            | **0.3455**                            |
| 2.5 | 0.2858            | 0.2413            | **0.2058**                            |
| 3.0 | 0.1562            | 0.1185            | **0.0922**                            |

**PED wins again**, once matched correctly -- the M=8 comparison in the
original notes was misleadingly easy to beat. So the desired ordering
(trace < static < PED) has now been checked and fails at BOTH ends of
the range tried (t=8 and t=40), not just t=40. This looks like a
structural property of this (128,64) SCL4 code/decoder configuration,
not something that improves by picking a different checkpoint within it.

**Where this leaves the drill**: within this code, checkpoint-based
pruning-from-64 (static or trace-ranked) has not beaten a
complexity-matched smaller full ensemble anywhere tested (t=8, 40, 48,
64). The trace-vs-static gain is real and reproducible at low-to-mid t,
but it's a gain against a weaker baseline (static_pm_rank), not against
the actually-relevant compute-matched alternative (PED). Candidate next
directions, not yet tried:
  - A different (n,k) code / SCL list_size -- Golay24 showed genuine,
    clean early-pruning gains over a full-PED-equivalent baseline earlier
    this project; this eBCH(128,64)/SCL4 configuration may simply not be
    a favorable case for the underlying idea.
  - Sweep the pruning ratio k/M itself (not just checkpoint t) at a few
    t values, in case a more aggressive or gentler prune ratio changes
    the picture -- untested so far, only k=8 has been tried throughout.
  - Reconsider whether SCL4's list diversity (only 4 competing paths per
    branch) is too narrow for AED-branch pruning to have real headroom
    over just running fewer, fully-decoded branches -- worth checking a
    wider list_size (SCL8/SCL16) where each branch's own decode is
    stronger and pruning might matter less, or the opposite (SCL2, where
    per-branch decode is weaker and ensemble diversity should matter more).

## Plan

1. Reuse `results/clean_trace_12k_upto_t64.csv` (already generated,
   cumulative through t64) and the unmodified
   `architectures/trace/train_mclass_trace_ranker.py` /
   `export_model_weights.py` -- no new dataset generation needed to
   start.
2. First experiment (cheapest, zero new code): retrain at t=64 with
   `--target log_gap` instead of `--target basin`, same hidden=64/epochs=20,
   compare against the same held-out split's basin-hit-rate metric the
   training script already reports, before spending compute on a full
   BLER sweep.
3. If log_gap alone doesn't close the gap, layer in: more epochs, a
   wider hidden layer, then finally multi-SNR training data.
4. Only promote a candidate to a real BLER sweep (via
   `common/mclass_bler_sim.cpp`, matching the exact t=64/k=8 slice above)
   once its held-out ranking metric clearly beats the current model's.
5. Keep the shipped `results/clean_trace_upto_t64_basin_h64_e20.pt`
   untouched throughout -- this drill produces new, separately-named
   models here, promoted to results/ only if they actually win.
