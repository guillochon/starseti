# Plan for responding to review 1

*Written 2026-09-24 against `paper/main.tex` (working copy, commit 7cb3be7 plus
the uncommitted edits the referee read). Each referee point was checked before
a fix was planned: the physics against a recomputation, the literature against
the Valency index, the statistics claims against `src/` and the `research/`
notes, the figure remarks against the PNGs. No edit has been made. Each item
below ends with a proposed action; approve, amend or strike them one at a
time.*

Verdict counts: 3 major, all legitimate. 9 medium: 7 legitimate, 2 legitimate
in substance but with an error in the referee's own arithmetic or reading that
the response should correct (m4, m5). 20 minor: 19 legitimate, 1 optional
(minor 14). Two of the referee's numbers are wrong (the RRab Kelvin--Helmholtz
time in m5; the reading of the duty cycle in m4), and one point (M1) has a
fresh observational anchor the referee did not mention (KM3NeT's 220 PeV
event, 2025) that softens its consequences for the RR Lyrae case.

A constraint on every abstract edit: the abstract is at 1919 characters
against the 1920-character, macro-free limit, so every change there must be
net-neutral or shorter.

---

## Major

### M1. Cross-section extrapolation — LEGITIMATE, load-bearing

**Check.** The paper's $E_{\tau=1}$ energies assume $\sigma\propto E$ from
1 TeV upward. Gandhi et al. (1996, Astropart. Phys. 5, 81, [hep-ph/9512364](https://arxiv.org/abs/hep-ph/9512364);
1998, PRD 58, 093009, [hep-ph/9807264](https://arxiv.org/abs/hep-ph/9807264))
give the deep-inelastic cross-section turning over to
$\sigma_{\rm CC+NC}\simeq7.8\times10^{-36}\,(E/{\rm GeV})^{0.363}$ cm$^2$
above $\sim$10 PeV, with a gradual softening from a few TeV. Both papers we
cite for "rises with energy" (Bustamante & Connolly 2019, [1711.11043](https://arxiv.org/abs/1711.11043);
IceCube 2021, [2011.03560](https://arxiv.org/abs/2011.03560)) measure the
sublinear regime, 18 TeV--10 PeV. Recomputed with the power law:

| path | $\Sigma$ (g cm$^{-2}$) | $\sigma$ for $\tau=1$ (cm$^2$) | paper | power law |
|---|---:|---:|---:|---:|
| Cepheid, central | $2.6\times10^{8}$ | $6.4\times10^{-33}$ | 912 TeV | $\sim$100 PeV |
| RRab, central | $3.7\times10^{9}$ | $4.5\times10^{-34}$ | 65 TeV | $\sim$70 TeV (at the crossover) |
| $\beta$ Cep, central | $2.7\times10^{10}$ | $6\times10^{-35}$ | 9 TeV | 9 TeV (linear regime) |
| Cepheid, grazing | $6.5\times10^{5}$ | $2.6\times10^{-30}$ | 2.5 EeV | $>10^{15}$ GeV, unreachable |
| RRab, grazing | $1.8\times10^{7}$ | $9.3\times10^{-32}$ | 640 PeV | $\sim$$2\times10^{20}$ eV |
| $\beta$ Cep, grazing | $4.2\times10^{8}$ | $4\times10^{-33}$ | 0.57 PeV | $\sim$30 PeV |

Every one of the referee's entries reproduces. The referee's internal-consistency
point (item 4) is also right and is independent of the law:
$\Sigma_{\rm disk}/\Sigma_{\rm graze}=2/(\sqrt{2R/H}f_{\rm above})$ is 401,
205 and 63 for the three classes, so under *any* single monotonic $\sigma(E)$
the central 912 TeV / 65 TeV / 9 TeV map to 366 PeV / 13 PeV / 570 TeV. Only
the $\beta$ Cep row of Table 1 matches; the Cepheid and RRab annular energies
in `research/modulation-physics.md` §5 were not computed from the same law
as the central ones. That is an error in the note, carried into the table and
`src/physics.py` (`beam_tau1_annular`).

**Consequences at the highest energy that is physical.** At 1 EeV,
$\sigma=1.4\times10^{-32}$ cm$^2$:

| class, annular | $f_{\rm abs}$ | $P_{\rm beam}/L$ ($\dtp=10^{-3}$, $d=0.5$) | gain over central |
|---|---:|---:|---:|
| Cepheid | $5.6\times10^{-3}$ | 0.09 | $\sim$15$\times$ |
| RRab | 0.16 | $3.2\times10^{-3}$ | $\sim$260$\times$ |
| $\beta$ Cep (at 30 PeV, $\tau=1$) | 0.6 | $8.3\times10^{-4}$ | $10^{3}\times$ |

Two anchors the referee did not use and the response should:

- **Glashow resonance.** A monochromatic $\bar\nu_e$ beam at 6.3 PeV has
  $\sigma_{\bar\nu_e e}\simeq5\times10^{-31}$ cm$^2$ per electron; in an
  H/He envelope ($n_e\simeq0.85/m_p$) that is $\Lambda\simeq3.9\times10^{6}$
  g cm$^{-2}$. On the grazing Cepheid column $\tau=0.17$, so
  $f_{\rm abs}\simeq0.15$ and $P_{\rm beam}/L\simeq3\times10^{-3}$: the
  Cepheid annulus is rescued, at a factor 4 rather than 100 from the old
  number. On the RRab ($\tau\simeq5$) and $\beta$ Cep ($\tau\simeq100$)
  grazing columns the resonant beam is absorbed *above* the zone, so it does
  not help those classes and they keep the tunable-energy design. IceCube
  detected one such event ([2110.15051](https://arxiv.org/abs/2110.15051),
  Nature 591, 220, 2021). Cost: two sentences plus a citation, as the referee
  says.
- **KM3NeT.** "All energies that nature demonstrably produces" can now be
  anchored at $\sim$220 PeV (KM3NeT Collaboration 2025, Nature 638, 376,
  [doi:10.1038/s41586-024-08543-1](https://doi.org/10.1038/s41586-024-08543-1)).
  At 220 PeV, $\sigma\simeq8\times10^{-33}$ cm$^2$: RRab grazing
  $f_{\rm abs}\simeq0.09$, Cepheid grazing $3\times10^{-3}$. So the
  $\beta$ Cep annulus (30 PeV) sits inside the observed spectrum, the RRab
  annulus is within a factor of a few of it, and the Cepheid annulus is not.

**Proposed actions (M1-a to M1-f):**

- **M1-a.** Section 2.4: quote the Gandhi et al. parameterization
  (linear to a few TeV, $E^{0.363}$ above $\sim$10 PeV), cite both papers,
  and recompute the three central $E_{\tau=1}$ (100 PeV, 70 TeV, 9 TeV).
  Replace "all energies that nature demonstrably produces" with the KM3NeT
  anchor and note that the Cepheid central beam sits at its edge.
- **M1-b.** Section 2.5 and Table 1: recompute the annular column from the
  consistent ratio, give the annular $E_{\tau=1}$ as 30 PeV for $\beta$ Cep
  and as "unreachable at $\tau=1$" for the other two, and add per-class
  $f_{\rm abs}$ and $P_{\rm beam}/L$ columns evaluated at the highest
  physical energy (propose 1 EeV, stated as the convention, with the 220 PeV
  KM3NeT value in the text). Add the Glashow line for the Cepheid. Drop
  "for every class" and "$8\times10^{-4}$ for every class".
- **M1-c.** Section 2.6: the Dyson-sphere fractions become per class
  (RRab: 3 per cent at $\eta=0.1$, 30 per cent at 0.01; Cepheid: most of a
  sphere at $\eta=0.1$ via 1 EeV, or 3 per cent via Glashow; $\beta$ Cep
  unchanged). The central-beam impossibility statement survives unchanged,
  because the central $\tau=1$ energies are all physical.
- **M1-d.** Section 2.7: the class ranking sharpens rather than weakens
  ($\beta$ Cep and RRab in the transmitter's favour, Cepheids against). One
  sentence.
- **M1-e.** Kick ceiling (Section 2.6) and Section 6.4: recompute with the
  RRab $f_{\rm abs}f_{\rm zone}$ at the chosen energy. At 1 EeV ($0.16$) the
  RRab energies per bit rise by $3.75\times$, $P_{\rm beam}/L$ to
  $(0.75$--$1.2)\times10^{-2}$, the Dyson share to 8--12 per cent at
  $\eta=0.1$, and the kick ceiling falls to $8\times10^{-5}$, so the
  measured thresholds are 30--50 times above it rather than 8--13. This is
  the largest downstream change and it reaches the abstract, Conclusion 3,
  Figure 3's dash-dotted line (`make_completeness.py`) and the viewer's
  physics box (`src/physics.py`: `F_ABS_TAU1`, `beam_tau1_*`, and the
  coupling must become per class). Decision needed: quote the RRab numbers
  at 1 EeV (conservative) or at 220 PeV (demonstrated); I recommend 1 EeV
  in the tables and 220 PeV in a clause.
- **M1-f.** Abstract: "1000 times less power than a central beam" becomes
  a per-class statement; shortest form "needs between 15 and 1000 times less
  power than a central beam, depending on class" (fits the length budget
  only if M2/m2's abstract edits are net-neutral; see below). Conclusion 2
  rewritten to match.
- **M1-g.** `research/modulation-physics.md`: append a dated correction
  section rather than rewriting history, since the paper says the derivation
  is in the released notes.

### M2. Sensitivity unmeasured for the ground-based 82 per cent — LEGITIMATE

**Check.** `research/injection-sensitivity.md` ("What this does not cover")
confirms it: no sparse injection has been run; the ASAS-SN and ZTF caches
now exist (`lc_cache.py`, `cache_index`), so the grid is runnable, and the
same note records that the sparse path needed a 10$\times$ larger injection
fraction to see anything. Section 6.4 and Conclusion 3 currently read as
statements about the whole search.

**Proposed actions:**

- **M2-a (preferred, needs compute).** Run `inject --find-threshold` on a
  few tens of ASAS-SN and ZTF RR Lyrae hosts (`--survey asassn` /
  `--survey ztf`, after the prefetch pass), staircase waveform, the 8-bit
  random word and the 28-bit primes word, at one bit per cycle *and* at the
  slow bit periods the epoch fold now searches (a code held for a night or a
  week per bit is the one the block path can read;
  `research/kick-and-epoch-routes.md`). Add the rows to Table 3 and a fourth
  panel or a second row to Figure 3. Expect bounds, not thresholds, at one
  bit per cycle; the useful number is the slow-code threshold.
- **M2-b (fallback, text only).** If M2-a cannot be run before resubmission,
  state in the abstract ("sensitivity measured on the TESS arm"), Section
  4.12, Section 6.4 and Conclusion 3 that the completeness is measured only
  for per-cycle TESS timing and the ground-based arm's contribution to the
  null is unquantified. The abstract already says "Injections into real TESS
  light curves", so the abstract cost can be a few characters.
- **M2-c.** Either way, Section 3.1's Gaia DR4 paragraph should say the
  block-averaged path's sensitivity is (measured / unmeasured), since it is
  what the DR4 extension rests on.

### M3. The zero-false-promotion control is a reproducibility check — LEGITIMATE

**Check.** Section 4.12 says the uninjected re-runs give "the false-alarm
rate"; the injection note says explicitly they are "a reproducibility check
and not a false-positive rate", and the hosts are drawn from stage-1 quiet
stars. The referee's replacement statistic is already in Table 2; per class
and survey the stage-1 flag rate is

| survey | Cepheids | RR Lyrae | $\beta$ Cep |
|---|---:|---:|---:|
| TESS | 7% | 21% | 55% |
| ASAS-SN | 14% | 46% | 12% |
| ZTF | 24% | 44% | (0 of 8) |

against gates set at $p<10^{-3}$, i.e. the white permutation null is wrong
for most real stars, which Section 3.3 already predicts.

**Proposed actions:**

- **M3-a.** Section 4.12: replace "which gives the false-alarm rate" with
  "which confirms that the quiet verdicts are reproducible under the
  injection machinery"; same fix in Section 4.9 ("at a measured cost of zero
  false promotions on 80 uninjected controls" becomes "with no quiet control
  changing verdict").
- **M3-b.** Add the table above (or one sentence with its range) to
  Section 4.9 or 6 as the empirical stage-1 false-positive rate, and say it
  is what the 11,102 SURVIVORs must be read against.
- **M3-c.** One sentence after Figure 3: completeness is conditional on the
  host being quiet; a code in a Blazhko star or a beating $\beta$ Cep is not
  represented. (`--any-host` exists and could measure this; propose as
  future work unless M2-a is run, in which case a small `--any-host` cell is
  cheap to add.)

---

## Medium

### m1. One-sided trigger asserted, not derived; pipeline is polarity-free — LEGITIMATE

**Check.** Section 2.1 ends "a sign convention that a search can exploit."
`research/sequence-pvalue.md` §Polarity (2026-09-22) records the opposite
decision: the matched filter tries both polarities and charges $\times2$;
seven of ten words in Table 5 are "inverted"; the decode is polarity-free.
The referee's physics is sound: a phase response curve of a nonlinear
limit-cycle oscillator is generically biphasic, and no Cepheid PRC has been
computed (Valency finds nothing on it). Option (b), requiring the sign, would
halve Equation 8's trials and remove the factor 2 in $T$, but the pipeline's
own note gives the reason it was not done: which symbol the transmitter kicks
on is its choice.

**Proposed action.** Option (a): hedge Section 2.1 to "we adopt the
one-sided model of Learned et al.; whether the phase response of a
$\kappa$-driven limit cycle is genuinely one-signed has not been computed,
and the search does not depend on the sign", delete "a sign convention that
a search can exploit", and add one sentence in Section 4.7 saying the
matched filter tries both polarities for that reason (it already charges for
it). Conclusion 1 keeps "one-sided" as the adopted model, with "as
proposed" added.

### m2. Abstract overstates SURVIVOR — LEGITIMATE

**Check.** Section 4.10 gives ground-based rows "the subset that applies";
Section 6 says most survivors reach the class "on excess per-epoch jitter
alone". The abstract's "survive every artifact, red-noise, systematics and
astrophysical veto" is not what those rows received.

**Proposed action.** Use the referee's wording, trimmed to fit the length
budget: "11,102 carry timing structure the vetting could not attribute to a
known cause, 97 per cent of them on block-averaged ground-based timing to
which only part of the vetting applies." (Net +30 characters; must be paid
for elsewhere in the abstract, e.g. by shortening the Learned sentence.)

### m3. The FDR pass cannot pass the injected positive control — LEGITIMATE, important

**Check.** The synthetic primes control names at $p=2.0\times10^{-3}$ (Table
5 note; `sequence-pvalue.md` §Polarity). The coincidence term for an exact
28-bit numeric hit is $2T\,R/\binom{L}{w}$ with $T\sim$ thousands, so no
exact 28-bit catalog word can reach $5\times10^{-7}$; roughly 40 bits would
be needed, which is the top of the $12\le L\le40$ window. As constructed,
the across-star FDR pass is therefore a test that a Learned-scale beacon of
the injected kind would fail, and the paper does not say so.

**Proposed actions:**

- **m3-a.** Say it plainly in Section 4.8: the FDR pass is a statement about
  *named-word claims*, a real 28-bit beacon would not pass it, and detection
  rests on the per-star gate plus the stage-2 requirements (the decoder finds
  the word on its own, squareness, the same word in each half of the
  baseline) that Section 6.3 lists. Add the minimum passable length.
- **m3-b (optional).** Note that a prior weighted toward prefixes of the
  short catalogs would let a primes beacon pass while changing little for the
  72--76 family; propose it as future work rather than changing the prior
  post hoc, since the store has been searched under the current one.

### m4. Kick ceiling duty cycle — LEGITIMATE as a clarity problem, but the referee's reading is not the only one

**Check.** Equation 5 with $d=0.5$ gives $6\times10^{-4}$; the paper's
$3\times10^{-4}$ is the $d=1$ value, and `physics.kick_ceiling` sets
`duty=1.0` deliberately: with no energy store the constraint is *per kick*
(one cycle's intake pays for one kick, whatever the average duty), so the
sentence "kicking every cycle with what it collected in that cycle" is
self-consistent. The referee is right that this is not explained and that
the same equation is used with $d=0.5$ everywhere else.

**Proposed action.** Keep $3\times10^{-4}$ (before the M1-e correction) but
rewrite the sentence: the kick ceiling is a per-kick limit, $\eta f_{\rm
coll}Gf_{\rm abs}f_{\rm zone}$ with no $d$, and a transmitter able to bank
even one cycle's intake doubles it for a fair-coin word. Then apply M1-e's
per-class $f_{\rm abs}f_{\rm zone}$. State in the response why $d=1$ is the
right reading. Alternative if the author prefers: adopt $6\times10^{-4}$
throughout, in which case Sections 6.4 and Conclusion 3 change as the
referee says.

### m5. Kelvin--Helmholtz times — LEGITIMATE mislabel; the referee's RRab number is wrong

**Check.** $GM^2/(RL)$ for the reference stars is $5.2\times10^{3}$ yr
(Cepheid) and $6.2\times10^{4}$ yr (RRab), not the referee's
$8\times10^{3}$ yr for the RRab. The paper's "1700 yr" and
"$2\times10^{4}$ yr" are exactly $\tau_{\rm KH}/3$ for both stars: the text
quotes the cancellation time $\tau_{\rm KH}/3$ and labels it $\tau_{\rm KH}$.
No envelope mass fraction was used; the referee's "goes the wrong way" rests
on their own arithmetic error.

**Proposed action.** Section 2.3: "$\tau_{\rm KH}=GM^2/(RL)$, about
$5\times10^{3}$ yr for the Cepheid and $6\times10^{4}$ yr for the RRab, so
the cancellation time $\tau_{\rm KH}/3$ is $\sim$1700 yr and
$\sim$$2\times10^{4}$ yr." Correct the referee politely in the response.

### m6. Squareness 0.93 versus 0.87 — LEGITIMATE, label them

**Check.** `research/irrational-beating-trap.md` table: 0.93 is the
injected counting-word staircase in the *synthetic RRab* (106 cycles); 0.87
is the primes control in the *synthetic Cepheid* (Sections 6.2, 6.3). Both
are synthetic controls of Figure 4, different stars and words.

**Proposed action.** Section 4.11: "the injected staircase code of the
synthetic RRab keeps its squareness at 0.93 (the synthetic Cepheid's primes
word, 0.87)".

### m7. Weighted Bonferroni is not Holm — LEGITIMATE

**Check.** The procedure in Section 4.8 is single-step with weights summing
to one; that is the union bound. Holm (1979) does contain a weighted
sequentially rejective variant, but the paper's procedure is not step-down.

**Proposed action.** Cite the union bound directly ("weighted Bonferroni"),
keep Holm only if pointing to the weighted section, and note in the response
that step-down Holm would gain little here because a single hit is being
tested per star. (`research/sequence-pvalue.md` uses the same citation and
should carry the same fix.)

### m8. Two missing anchors — LEGITIMATE

**Check.** (a) Edmondson & Stevens (2003, IJA 2, 231, "The utilization of
pulsars as SETI beacons") is not on arXiv, so Valency cannot confirm it;
ADS must. The lineage claim is right: Edmondson is a Hippke et al. (2015)
co-author. While checking, found that `references.bib` spells the fifth
Hippke et al. author "Steven, Ian R."; it is Stevens. (b) Gandhi et al. 1996
and 1998 confirmed via Valency (DOIs 10.1016/0927-6505(96)00008-4 and
10.1103/PhysRevD.58.093009).

**Proposed actions.** Add both Gandhi entries and Edmondson & Stevens 2003
to `references.bib` (Edmondson via ADS), cite E&S in Section 1.1 before
Chennamangalam et al., cite Gandhi in Section 2.4 (part of M1-a), and fix
the Stevens typo.

### m9. Which series the whiteness statistics read — LEGITIMATE

**Check.** `analyse_oc` (`stage1_pulse.py` ~3746--3800): the ACF,
Ljung--Box and permutation entropy run on the detrended *level* series
`oc_d`; the decode statistics run on `code_increments` and the shuffle null
for those is differenced per draw. Under the staircase model a fair-coin code
makes the level series a random walk, which the level Ljung--Box flags (this
is the excursion that carried the 28-bit primes word to triage at
$2.4\times10^{-3}$ in the injection note), while the step series is white.

**Proposed action.** One sentence in Section 4.6 ("computed on the detrended
level series") and one in 4.7 ("all four read the steps"), plus a clause
saying the two are complementary: the level statistics see a code's
excursion, the step statistics see its bits. The injection note's
"concrete next step" (a Ljung--Box on the increments) can be named as future
work.

---

## Minor

All checked; actions are one clause each unless noted.

1. Abstract "1000 times" — follows M1-f.
2. "the authors" — appears in Section 5 ("the viewer the authors use") and
   the acknowledgments. Change Section 5 to "the author". The AI-disclosure
   paragraph is fixed text and stays as is; say so in the response.
3. $\Delta P/P\simeq3\Delta E/(GM^2/R)$ — verified: $P\propto R^{3/2}$ gives
   $\Delta P/P=\tfrac32\Delta R/R$, and the virial energy $E=-GM^2/(2R)$
   gives $\Delta R/R=2\Delta E/(GM^2/R)$. Add the clause.
4. "subtends $2.7^\circ$" is the angular diameter (the note's table says
   2.66$^\circ$, half-angle 1.33$^\circ$). Say "angular diameter".
5. Rimoldini et al. classify 12.4 million *candidates*. Change the word.
6. Cepheid window 0.05--60 d (`StarClass("cepheid", 0.05, 60.0, ...)`).
   Anomalous Cepheids have $P\gtrsim0.3$ d, so 0.05 d does admit $\delta$ Sct
   misclassifications. Add a clause: the Gaia period is used within
   $\pm2$ per cent where it exists, so the wide window binds only for the
   few rows without one; or narrow it (code change, then a clause).
   Recommend the clause.
7. TESS $\beta$ Cep per-cycle precision: 129 s median formal error,
   822 s O--C rms (injection note). Add to Section 4.4.
8. Small-$\delta t$ linearization in Equation 7: add a clause; the steep
   rising branch of an RRab is where $\dot F$ changes fastest and where a
   large $\delta t$ biases the implied offset low. Legitimate.
9. Otsu $2/\pi$: add "(the between-class share for a normal split at its
   mean, $2\phi(0)^2/(\tfrac12\cdot\tfrac12)=2/\pi$)". Legitimate.
10. The 1 per cent for a 21-bit fair-coin rotation word was by simulation:
    20,000 draws (0.9 per cent at 21 bits; `irrational-beating-trap.md`).
    Say "by simulation".
11. Table 2 Kepler rows: footnote that the Kepler RR Lyrae field sample is
    not included and is an obvious dense-mode extension.
12. Table 4 "BCEP?": footnote that the query mark is Gaia's and that the
    text argues the star is SPB-like.
13. Section 6.2: 93 of 1123 named words is 8.3 per cent; add it.
14. Mass functions: optional rewording; recommend keeping the numbers, since
    the referee checked them and they are the evidence.
15. Section 7: add that the timing channel is one-sided and peak-limited,
    for which the capacity-achieving input is discrete rather than Gaussian
    (Smith 1971, Inf. Control 18, 203; not on arXiv, verify via ADS), so a
    transmitter restricted to fixed kicks is *closer* to the two-level case
    than the Gaussian argument suggests; the qualitative conclusion stands.
16. Figure 1(a): `make_geometry.py` draws the zone at $r=0.86R$ ("exaggerated
    inward" in the script); add "zone depth exaggerated" to the caption.
17. Figure 3: legend says "triage", caption says "flagged" (confirmed in
    `make_completeness.py` line 80); use "flagged". The "Learned et al.
    $10^{-3}$" label overlaps the RR Lyrae curves (confirmed in the PNG);
    move it to the top of the panel. The dash-dotted kick-ceiling line moves
    with M1-e.
18. Figure 4(d) right: drop the named window's $R^2$ from the legend
    (`make_candidates.py` line 299) and keep it in the text.
19. References: `paper/README.md` says volume and page fields from Valency
    are unverified while the acknowledgments say every entry was validated
    against ADS. Do the validation before resubmission (this is also where
    Edmondson & Stevens and Smith 1971 come from), then keep the sentence
    truthfully.
20. "per cent" (about 25 uses) versus "%" (about 16 uses with digits):
    pick "per cent" in prose and "%" in tables, and apply throughout.

---

## Order of work

1. **Physics block first** (M1-a..g, m4, m5, m8b), since it changes numbers
   that appear in the abstract, Table 1, Sections 2.4--2.7 and 6.4,
   Conclusions 2--3, Figure 3 and `src/physics.py`, and every later edit
   should be made against the corrected numbers. Decide the "highest
   physical energy" convention (1 EeV vs 220 PeV) before starting.
2. **Sparse injection grid** (M2-a) started in parallel, because it is the
   only item with wall-clock cost; text falls back to M2-b if it does not
   finish.
3. **Statistics wording** (M3, m1, m3, m7, m9), each a paragraph or less.
4. **Minors and figures** (m2, m6, minors 1--20), then the ADS pass over the
   bibliography (minor 19) last, so new entries are covered.
5. **Response letter** drafted from this file, including the two corrections
   to the referee (m4 reading, m5 arithmetic).

## Files each item touches

| item | files |
|---|---|
| M1 | `paper/main.tex` (abstract, §2.4--2.7, Table 1, §6.4, Conclusions), `paper/references.bib`, `src/physics.py`, `paper/figures/make_completeness.py`, `research/modulation-physics.md` |
| M2 | injection run (`data/injection.sqlite`), `main.tex` (abstract, §3.1, §4.12, Table 3, §6.4, Conclusion 3), `make_completeness.py` |
| M3 | `main.tex` §4.9, §4.12, §6 |
| m1 | `main.tex` §2.1, §4.7, Conclusion 1 |
| m2 | `main.tex` abstract |
| m3 | `main.tex` §4.8 |
| m4 | `main.tex` §2.6 (and §6.4 via M1-e) |
| m5 | `main.tex` §2.3 |
| m6 | `main.tex` §4.11 |
| m7 | `main.tex` §4.8, `research/sequence-pvalue.md` |
| m8 | `main.tex` §1.1, §2.4, `references.bib` |
| m9 | `main.tex` §4.6, §4.7 |
| minors | `main.tex` throughout; `make_geometry.py`, `make_completeness.py`, `make_candidates.py` captions/legends; `references.bib` |

---

## Code runs to perform before the edits

*Added 2026-09-24, after checking `src/inject.py`'s parser, the injection
store and the cache.* Two facts change the plan above slightly:

- **A sparse ASAS-SN grid already exists**: run `V-asassn-stage2`
  (`data/results/injection/efficiency-V-asassn-stage2.md`), 40 hosts per
  class, staircase, `random:8` and `primes`, 1 and 2 cycles per bit, fixed
  $\dtp$ = 1e-3..3e-2, with stage 2. ASAS-SN RR Lyrae at one bit per cycle
  never reach 50 per cent triage below 3e-2 (at two cycles per bit the primes
  word reaches it at 2.3e-2). It was run between 2026-09-17 and 09-21 and so
  predates the kick and epoch-fold routes (commit 0f97f7e, 09-22), which are
  exactly the routes built to read a code on a block-averaged star. It is a
  usable *bound* today and the reason to re-run with the current code.
- **No prefetch is needed.** The sparse cache holds ~35k ASAS-SN and ~127k
  ZTF files, and `hosts_from_store` only draws cached quiet hosts; the V run
  proves ASAS-SN hosts are drawable. If a ZTF run reports too few hosts,
  prefetch with `stage1 --survey ztf --lc-cache-prefetch --targets ...`.

Each run reads `data/cepheid.sqlite` read-only and writes only
`data/injection.sqlite` and `data/results/injection/`, so they can run while
nothing else is writing the store. Order by value:

1. **M2-a, ASAS-SN RR Lyrae, current code, slow codes included.**
   `./starseti inject --survey asassn --classes rrlyrae --hosts 40 --words random:8,primes --bit-cycles 1,8,32 --find-threshold --run-id R1-asassn-rrlyrae --jobs 4`
   About 2,000 injections plus 40 controls (the V run was 1,920). If it is
   slow, drop to `--bit-cycles 1,8` and `--hosts 30`.
2. **M2-a, ZTF RR Lyrae, same grid.**
   `./starseti inject --survey ztf --classes rrlyrae --hosts 40 --words random:8,primes --bit-cycles 1,8,32 --find-threshold --run-id R2-ztf-rrlyrae --jobs 4`
3. **M3-c, completeness and false-promotion rate on the population as it is**
   (not only quiet hosts), TESS RR Lyrae at the two published thresholds:
   `./starseti inject --survey tess --classes rrlyrae --hosts 40 --any-host --words random:8,primes --frac 2.4e-3,4e-3 --run-id R3-tess-rrlyrae-anyhost --jobs 4`
   The uninjected control rows of this run are the first real measurement of
   how often a non-quiet host is promoted by the injection path, and the
   injected rows give completeness unconditional on quietness.
4. **Table 3 note, the two-mode beat control the paper says has not been
   run** (not asked by the referee, but promised in the manuscript):
   `./starseti inject --survey tess --classes rrlyrae --hosts 40 --words beat,beat:0.707 --snr 3,5,10 --run-id R4-tess-rrlyrae-beat --jobs 4`
5. **Optional, M2-a for the other sparse classes** (Cepheids are 20 per
   cent of the ASAS-SN arm): repeat run 1 with `--classes cepheid` and a
   new `--run-id`.
6. **After review of each report, publish the thresholds to the store** so
   the viewer's physics box stops falling back to the Learned scaling for
   the sparse rows:
   `./starseti inject --report-only --run-id R1-asassn-rrlyrae --publish`
   (and the same for R2). Do this only after the M1-e change to
   `src/physics.py`, or the box will convert the new thresholds with the old
   per-class coupling.

Runs 1--2 feed Table 3 and Figure 3; run 3 feeds M3-b/M3-c; run 4 feeds the
Table 3 note and Section 4.11. Nothing else in the plan needs compute before
the edits: the flag rates for M3-b come from Table 2, the physics numbers
for M1 are analytic, and the bibliography pass is ADS work.

## Status of the code side of M1 (2026-09-24)

Done, ahead of the paper edits, so that sparse thresholds can be published
with the corrected coupling:

- `src/physics.py`: cross-section model (`sigma_nu_n`, linear to the Gandhi
  et al. 1998 power law), columns and $\tau=1$ energies computed from the
  reference stars instead of stored as strings, per-class `f_abs` at
  `BEAM_ENERGY_MAX_GEV` = 1 EeV (with the KM3NeT 220 PeV value carried as
  `BEAM_ENERGY_OBSERVED_GEV`), a `glashow` geometry, `kelvin_helmholtz_yr`,
  and per-class kick ceilings in the report and ceiling tables. Couplings
  are now 1.67 ($\beta$ Cep), 6.5 (RRab), 178 (Cepheid) for the annulus;
  6.0 for the Glashow beam on the Cepheid; central unchanged.
- `src/candidates.html`: the physics-box tooltip states $f_{\rm abs}$ and the
  beam energy per class, and offers the Glashow line where it applies.
- `paper/figures/make_completeness.py`: one kick-ceiling line per panel (the
  PDF is not regenerated until the figure edits of minor 17 are approved).
- `README.md` physics paragraph; `tests/test_physics.py`, `test_routes.py`,
  `test_serve_candidates.py` updated to the per-class numbers.
- Not done: `research/modulation-physics.md` (M1-g) and every number in
  `paper/main.tex`, which wait for approval.

Consequence for the viewer: a Cepheid's energy-per-bit box rises by 100x
and an RR Lyrae's by 4x once the server restarts on this code.

## Status after the second pass (2026-09-24, evening)

The four injection runs came back and every item of the plan is implemented
except the two that need resources this session does not have. In brief:

- **Runs.** R1/R2 (sparse RR Lyrae): a one-bit-per-cycle code is invisible on
  ASAS-SN and ZTF (flag rate <= 18% below the 0.05 ceiling); a code at 32
  cycles per bit is flagged at 50% for dt/P = 6e-4 to 2.6e-3 and at 8 cycles
  per bit for 2e-3 to 6e-3; nothing sparse was ever decoded or named. R3
  (`--any-host`): 60% of the cached TESS RR Lyrae pool is flagged uninjected;
  decode completeness on those hosts is 65% (8-bit) and 31-41% (primes) at
  2.4e-3, i.e. no worse than on quiet hosts. R4 (beat control): 74-94% flagged,
  0 of 228 promoted by a decode route, 0 decoded, 8% named at some p, none
  below the gate. Details in `research/injection-sensitivity.md`.
- **Paper.** Every M, m and minor item is in `paper/main.tex`: the abstract
  (1912 characters), Sections 1.1, 2.1, 2.3-2.7, Table 1 (regrouped so it fits
  the text width), 3.1, 4.3-4.8, 4.9, 4.11, 4.12 (two new paragraphs), Table 3
  (18 rows, with the sparse arm), Figure 3 (two rows; legend says "flagged";
  the Learned label moved), Figure 4 (named-window R^2 dropped from the
  legend), Figure 1 caption, Table 2 and 4 notes, Sections 5, 6.2, 6.4, 7,
  Conclusions 1-3, and "per cent" throughout. Both figures were regenerated
  and the PDF builds cleanly (28 pages; the one stuck-float warning predates
  this revision).
- **Bibliography.** Gandhi 1996 and 1998, IceCube 2021 (Glashow), KM3NeT 2025,
  Smith 1971 and Edmondson & Stevens 2003 added; "Steven" fixed to "Stevens"
  in Hippke et al. (2015).
- **Notes.** `modulation-physics.md` carries a dated correction section (M1-g),
  `injection-sensitivity.md` the new runs, `sequence-pvalue.md` the Holm fix.
- **Open.** (1) Minor 19, the ADS pass over every bibliography entry, needs an
  ADS token or a manual session; the acknowledgment's sentence should not be
  left as is until it is done. (2) The sparse-arm thresholds are quoted from
  logistic fits that exclude the ceiling row; `inject --find-threshold`
  should do that itself before the numbers are published to the store with
  `--publish`, which has therefore not been run for R1/R2. (3) The response
  letter to the referee has not been drafted; this file is its source.
