# A Benjamini-Hochberg pass over the sequence matches

> Plan, not yet implemented. The open item at the end of
> [sequence-pvalue.md](sequence-pvalue.md) ("There is still no across-star
> correction") is what this note designs. Store numbers measured 2026-09-19.

`seq_p` is a look-elsewhere correction over the catalog search *within one
star*: every packed stream, rotation, offset and length the matcher tries is
already inside it. It says nothing about the fact that we ran that search
55,601 times. The 0.01 gate is therefore a per-star screening rule wearing the
clothes of a global statement, and since `is_followup` promotes on it, the
promoted set has no false-discovery rate attached to it at all.

## What the pass would say today: nothing survives

| `seq_p` bin | rows |
| --- | ---: |
| < 1e-3 | 1 |
| 1e-3 .. 1e-2 | 3 |
| 1e-2 .. 0.05 | 22 |
| 0.05 .. 0.1 | 29 |
| 0.1 .. 0.5 | 203 |
| 0.5 .. 1 | 217 |
| = 1 (unmatched) | 55,126 |

Four rows clear the 0.01 gate and **one of them is `SYNTH cepheid`**, an
injection control living in `stage1_sequence` next to the real stars. The three
real ones:

| star | label | L | `seq_p` | `seq_squareness` |
| --- | --- | ---: | ---: | ---: |
| ASASSN 575526037701 | counting [u, 43, 44, 45, 46, u] | 29 | 3.4e-4 | NULL |
| ASASSN 214749425367 | counting [u, 27, 28, 29, 30, u] | 39 | 3.6e-3 | NULL |
| ASASSN 335007454137 | counting [u, 6, 7, 8, 9, 10, u] | 26 | 5.8e-3 | NULL |

BH at q = 0.05 over the 41,844 real timed stars needs the smallest p below
0.05/m = 1.2e-6. The best real star is 3.4e-4, short by a factor of ~285. The
null result is the point: it is the statement the paper should be making, and
it is stronger than the uncorrected one it makes now.

Worth noticing separately that all three have a null `seq_squareness`, so they
reached the gate without the two-level-fit check that keeps smooth astrophysics
out. That is a question about the gate, not about the FDR pass, but it should
be answered before any of them appears in a candidate table.

> Answered (2026-09-19): the squareness was null because the fold's own word
> needs three repeats and these sit at 2.5. The catalog route now requires the
> *named window* to beat a sinusoid ([sequence-pvalue.md § The window has to
> be square](sequence-pvalue.md#the-window-has-to-be-square)); one of the three
> (575526037701, 0.30) survives it, the other two (0.00, 0.17) do not.

## The decisions

### 1. What counts as a test

Three exclusions before anything is counted:

- **`SYNTH` rows** (3 in `stage1_sequence`). Injections are not discoveries.
- **Errored stage-1 rows.** 55,601 sequence rows join to 41,847 rows with
  `error IS NULL`; only the latter were ever really tested.
- **Cross-survey duplicates.** One physical star arrives as up to three survey
  targets. Among the 472 matched rows there are 470 distinct `source_id`s, so
  the duplication is negligible at the top of the list but it inflates `m`.
  Deduplicate by `source_id`, take the minimum p per star, charge a within-star
  Bonferroni by the number of surveys that tested it.

### 2. Discreteness, which is where the power is

`seq_p` has a point mass at exactly 1.0 covering 99.1% of rows, because most
stars have no word long enough to be nameable at all. `seq_bit_len` is stored
for unmatched rows too (0..64), so the *minimum achievable* p-value is
computable for every star without redoing any matching:

| word length | stars | can ever be named? |
| --- | ---: | --- |
| L < 12 | 28,985 | no, p is identically 1 |
| L 12-16 | 4,531 | pi/e prefix only (numeric needs 17) |
| L >= 17 | 8,328 | either corpus |

Tarone's reduction drops tests whose minimum achievable p exceeds the
threshold. This preserves control and shrinks `m` from 41,844 to ~12,859, a
3.3x gain, and it is the standard move for discrete tests (Gilbert's discrete
BH). It does not rescue the present data -- the bar moves only to 3.9e-6 -- but
it is free and it matters more as the store grows.

### 3. Dependence

Multi-survey rows of one star are positively dependent. BH controls FDR under
positive regression dependency, which is defensible here. Benjamini-Yekutieli
is the conservative alternative and costs a factor of log(m) ~ 9.5, which is
not worth paying. Deduplicate first (decision 1), use BH, state the assumption.

### 4. The conservativeness trap

`seq_p` is *already* a weighted-Bonferroni FWER correction. Stacking BH on it is
valid but doubly conservative: under the null `seq_p` is stochastically much
larger than uniform, so Storey's pi0 comes out at ~1 and adaptive BH buys
nothing. The way to recover power is an **empirical null** -- the injection
machinery already runs quiet controls (320 of them, per
[injection-sensitivity.md](injection-sensitivity.md)), so the actual null
distribution of `seq_p` on real hosts with no code injected can be measured and
calibrated against, instead of the analytic bound. That is the difference
between "BH says no" and "BH says no, and here is what it could have seen".

Treat this as a separate follow-on: it needs a control-injection campaign and is
really a sensitivity measurement, not a correction.

### 5. Where it runs

`SEQ_P_GATE = 0.01` is consumed in three places: `score.sequence_bonus`
(`score.py:193`), the triage route in `is_followup` (`stage1_pulse.py:3566`),
and the viewer. A BH threshold is a population quantity, unknowable until every
star is timed -- and stage 1 evaluates `is_followup` per star while it is still
writing. **So BH cannot drive promotion.**

Keep 0.01 as the cheap per-star screening gate for triage; add BH as a
reporting layer over the finished store. There is precedent:
`store.recompute_triage` (`store.py:860`) already re-applies a predicate across
the whole store so that a gate change does not force a stage-1 rerun.

## Implementation sketch

- A `seq_q` column on `stage1_sequence` (the BH-adjusted q-value) **and** a
  `seq_q_m` recording the family size it was computed against. A q-value is
  meaningless without its `m`, and `m` changes the moment the faint ZTF sample
  lands.
- `store.recompute_sequence_fdr()` in the `recompute_triage` mould: read all
  `(target, survey, seq_p, seq_bit_len, source_id)`, drop SYNTH and errored rows,
  dedup by `source_id`, compute the minimum achievable p for the Tarone
  reduction, run BH, write back.
- Wire it as `starseti stage2 --fdr` or a small `stage3`. No network, no light
  curves: seconds, not hours.
- Viewer: show `seq_q` beside `seq_p` and make it filterable, so "survives FDR"
  is a visible property of a row rather than a claim in prose.
- Tests: a known p-vector against a reference BH implementation; the Tarone
  reduction preserving control on synthetic discrete nulls; the SYNTH-exclusion
  invariant.

## What it changes in the paper

- The named-word count in the sequence section is currently called a display
  quantity "until" this exists. It stops being one.
- "No star decodes to a recognizable sequence at p < 1e-3" becomes an FDR
  statement, which is both stronger and more defensible.
- It is a **precondition for the faint ZTF expansion.** Going to G <= 18.5 grows
  the ZTF RR Lyrae arm by up to 14x (135,297 Galactic RR Lyrae pass that cut,
  against 9,637 in the store today), which moves the BH bar by another order of
  magnitude. This wants to be in place before those rows land, not after.

## Order of work

Build decisions 1, 2, 3 and 5 first: unambiguously correct, roughly a day, and
it gives the paper a proper null. The empirical null (decision 4) is a separate
piece of work with its own compute cost.
