# Plan for responding to review 3

*Written 2026-09-26 against `paper/main.tex` (the working copy the referee
read, after commit c5569ca, plus the uncommitted text already in the tree).
Each point was checked against the source and, where a count was in
question, against `data/cepheid.sqlite` and `data/injection.sqlite`. No edit
to the paper has been made. Each item ends with a proposed action; approve,
amend or strike them one at a time.*

Verdict: a polish pass. Every consistency item is real. Two of the referee's
diagnoses are the wrong mechanism (1.3, and the accounting inside 1.7), and
one suggested sentence (the abstract) does not fit the 1920-character limit.
The live store still matches the paper's survivor split among currently
flagged rows. It is three flagged rows ahead of the paper (38,233 against
38,230), all three WEAK. That drift is not part of this pass.

The abstract is at 1920 of 1920. The one abstract edit below is shorter than
the sentence it replaces.

---

## 1. Consistency between sections

### 1.1 The 11,102 in Section 7 — LEGITIMATE

**Check.** Section 7 still says an efficient message could sit among the
11,102 SURVIVORs. The abstract, Section 6 and Conclusion 4 all say 6823.
11,102 is the first-pass count from before the ground re-vet.

**Proposed action.** Replace 11,102 with 6823 in Section 7.

### 1.2 The magnitude cut against the ZTF window — LEGITIMATE

**Check.** Section 4.1 applies \(G\le 15.5\) to the target list and then, two
sentences later, matches ZTF inside \(13\le G\le 18.5\). Section 3.1 repeats
the faint ZTF window. ASAS-SN has its own window too (\(10.5\le G\le 15.5\)).
The \(G\le 15.5\) cut is the TESS list. ZTF is matched to the catalog inside
the survey window, which is why it reaches 61,975.

**Proposed action.** In Section 4.1: "a magnitude cut of \(G\le 15.5\) is
applied for TESS; ZTF and ASAS-SN use their own windows below."

### 1.3 The Table 2 cadence split — LEGITIMATE, wrong mechanism

**Check.** The note's 4600 + 8293 and 1345 + 110 are one above the table
cells. They are not rows counted twice. The store's dense counts include the
three synthetic controls, one per class, all timed per cycle. The table
excludes them, as the note's first sentence says. Removing those three:

| class | two-minute, real stars | full-frame | sum | table |
|---|---:|---:|---:|---:|
| RR Lyrae | 4599 | 8293 | 12892 | 12892 |
| \(\beta\) Cephei | 1344 | 110 | 1454 | 1454 |
| Cepheid | 2702 | 57 | 2759 | 2759 |

The Cepheid split is not quoted in the note, and it already matches.

**Proposed action.** Quote 4599 and 8293, and 1344 and 110, so the split
sums to the table. One clause: the synthetic control of each class is in
the store's dense count and not in the table.

### 1.4 "Three further controls" — LEGITIMATE

**Check.** The last paragraph of Section 4.12 says three and then numbers
First through Fourth (quiet-host condition, survey flag rate, beat control,
thermal tail).

**Proposed action.** Change "Three" to "Four".

### 1.5 Uninjected controls, absolute then qualified — LEGITIMATE

**Check.** The opening of Section 4.12 says none of the 80–160 uninjected
controls per class reached the vetting gate. The thermal-tail paragraph of
the same section says six of that run's 80 did. Both are right about their
own runs. The first sentence reads as a claim about every run.

**Proposed action.** "none of the 80–160 uninjected controls per class in
the main runs reached the vetting gate".

### 1.6 The re-vet tally against the survivor tally — LEGITIMATE

**Check.** The four classes in the re-vet sentence sum to 15,796, not
15,798. The target list (`revet_ground_survivors.txt`) is no longer in the
tree, so the two missing rows have no class we can still name.

The 6,670 do reconcile with the present ground-based survivors, once the
rows that have left triage are put back:

- 5,979 ZTF survivors, all still flagged
- 639 ASAS-SN survivors still flagged
- 52 ASAS-SN survivors whose stage-1 row is no longer in triage (30 still
  carry a verdict, 22 do not)

5,979 + 639 + 52 = 6,670. The 197 TESS full-frame survivors were not in the
15,798. The present sparse total is \(6670 - 52 + 197 = 6815\), and with the
eight dense survivors that is 6823.

**Proposed action.** One clause after the re-vet sentence: of the 6,670,
52 ASAS-SN rows have since been re-timed out of triage, leaving the 5,979
ZTF and 639 ASAS-SN survivors still flagged; the 197 full-frame rows were
not in the 15,798. Write the class sum as 15,796 of the 15,798, and do not
invent a class for the other two.

The ZTF injection sentence is the same kind of gap, and the injection store
closes it. Run `R5-ztf-rrlyrae-stage2` has 414 flagged injections with a
stage-2 class. Seventeen of them sit at \(\delta t/P = 0.05\), the amplitude
at which the injection stops being a phase code. The 397 below that line are
86 SURVIVOR, 36 WEAK, 159 ARTIFACT, 113 ASTRO, 2 SYSTEMATICS and 1 RED NOISE.

**Proposed action.** Name the three: two SYSTEMATICS and one RED NOISE.

### 1.7 "Two of them" and the eight dense survivors — LEGITIMATE

**Check.** The eight dense SURVIVORs, in score order, are the six rows of
Table 5 plus two the text never names. The referee's parenthetical counts
nine, because BD+46 3337 is one of the three \(\beta\) Cephei stars in the
table, not an extra star beside three beats. The table's six are five with
an explanation and one without:

1. TIC 152825245, RRc, alias of a prewhitened 0.61 mode
2. TIC 240268877 (BD+46 3337), SPB-like
3. TIC 10931356, RRab, unexplained
4. TIC 467480241, \(\beta\) Cephei beat, 3.4 d
5. TIC 114228663, \(\beta\) Cephei beat, 3.2 d
6. TIC 142240683, RRc, alias of a prewhitened 0.61 mode

The other two, both annotated "unresolved beat" and both with
`common_frequency` false, so the hard veto never fires:

- TIC 235416406, \(\beta\) Cephei, \(P = 0.116\) d, 845 cycles, score 0.79.
  Five modes prewhitened at 0.8–3 per cent of the dominant amplitude. The
  O–C peaks at 0.032 cycles\(^{-1}\) and the amplitude at \(5\times 10^{-4}\).
  Secondary periods at 3.6 d and 3.1 d.
- TIC 347799445, RRc, \(P = 0.280\) d, 268 cycles, score 0.76. A 0.171 d
  mode (frequency ratio 1.64, the 0.61 family) was prewhitened. The O–C
  peaks at 0.220 cycles\(^{-1}\) and the amplitude at 0.0025. Its O–C
  scatter is 108 s against a median formal error of 851 s, the same
  inflation discussed under the sentences item for TIC 147576267.

"Two of them" has to mean these two. The alias route is already named, for
the two RRc stars in the table, and TIC 10931356 is already the unexplained
one. Conclusion 4's "aliases and mode beating, with one RRab unexplained"
covers the table and leaves these two out.

**Proposed action.** In the last paragraph of Section 6.1, name the two and
say why they survive: the common-frequency veto requires the O–C and the
amplitude to peak together, and in both stars they do not. Add the same
half-sentence to Conclusion 4.

### 1.8 "The same window" — LEGITIMATE

**Check.** Table 6's four counting matches are 73–76, 72–75, 73–75 and
72–75. They overlap on 73–75. They are not one window.

**Proposed action.** "four overlapping windows of the integers 72–76".

### 1.9 "Stage 0" — LEGITIMATE

**Check.** Section 3.1 says the classifications and periods are what stage 0
consumes. Section 4 introduces stages 1 and 2 and describes the target list
without giving it a number. The store table is `stage0_targets`.

**Proposed action.** In the opening paragraph of Section 4, call the
target-list build stage 0.

### 1.10 The 25 per cent and the 60 per cent — LEGITIMATE

**Check.** Table 2's TESS RR Lyrae flag rate is \(2684/12892 = 21\) per
cent, because the cell mixes two-minute and full-frame rows. Section 4.3
states the 25 per cent for the two-minute rows alone, against 54 per cent
for the rows the two per cent floor had timed. The 60 per cent in Section
4.12 is the separate pool of 40 hosts, timed with that old floor. The
parenthetical already says the pool was timed with the old floor, and then
cites Table 2 for the 25 per cent.

**Proposed action.** Cite Section 4.3 for the 25 per cent, and add that the
60 per cent is this pool of 40, not the 54 per cent of all first-pass
two-minute rows.

---

## 2. Notation collisions

### 2.1 \(G\) the gain and \(G\) the gravitational constant — LEGITIMATE, other symbol

**Check.** Newton's constant appears twice, both in Section 2.3
(\(E\simeq -GM^2/2R\) and \(\tau_{\rm KH}=GM^2/RL\)). It is not in the
Table 1 note; the \(G=1\) there is the gain. The gain is the paper's symbol
from Equation 2 on, including the phrase "\(G=1\) kick ceiling", which the
last revision just established in the text, the figure and `src/physics.py`.

**Proposed action.** Write Newton's constant as \(G_{\rm N}\) in Section 2.3.
Leave the gain as \(G\). Renaming the gain would touch every ceiling
statement for a collision that lives in one subsection. DECISION if you
would rather follow the referee and rename the gain; the recommendation is
\(G_{\rm N}\).

### 2.2 \(\eta\) the efficiency and \(\eta\) the growth rate — LEGITIMATE

**Check.** The growth rate is one sentence in Section 2.4
(\(W = 2\eta E_{\rm puls}\)). Everywhere else \(\eta\) is the collector
efficiency (Section 2.6, Equation 4, Figure 1).

**Proposed action.** Use \(\kappa\) for the growth rate in that sentence.

### 2.3 \(L\) the luminosity and \(L\) the word length — LEGITIMATE

**Check.** Luminosity \(L\) is Section 2. Word length \(L\) starts at
Section 4.6 (\(L\ln 2\), \(12\le L\le 40\)), is a column of Table 6, and is
in the Figure 4 caption. The two do not share an equation. \(\ell\) is
already the spherical-harmonic degree in Section 2.5, so it is the wrong
replacement.

**Proposed action.** Rename the word length to \(N_{\rm bit}\) in Section
4.6, the coincidence formula, Table 6 and the Figure 4 caption. Leave
luminosity as \(L\).

### 2.4 \(d\) the duty cycle and d the day — LEAVE

**Check.** The referee is right that a reader can stumble, and also right
that the collision is conventional. In the current text the duty cycle is
introduced as "duty cycle \(d\)", and "0.57 d" is the reference period in
the Section 2 preamble, not in the same sentence as \(d = 0.5\).

**Proposed action.** No edit.

### 2.5 \(T\) the TESS magnitude — LEGITIMATE

**Check.** Section 6.3 writes \(T = 12.45\), \(G = 12.60\) for
TIC 152825245, then \(T = 13.21\) for TIC 147576267, with no expansion.
\(G\) on that page is the Gaia magnitude, which Section 4.1 has already
established. \(T\) has not been.

**Proposed action.** At the first use: "TESS magnitude \(T = 12.45\), Gaia
\(G = 12.60\)".

---

## 3. Sentences and phrases

### 3.1 Abstract, the RR Lyrae star as the agent — LEGITIMATE, wording has to be shorter

**Check.** "an RR Lyrae star can drive it only with most of its light" makes
the star the transmitter. Section 2.6 says the beam needs most of the star's
light. The referee's replacement ("and around an RR Lyrae star it needs most
of the star's light to power") is 9 characters longer than the clause it
replaces, and the abstract is at 1920 of 1920.

**Proposed action.** "driving an RR Lyrae star takes most of its light",
which is 9 characters shorter (abstract 1911). Run `util/check_abstract.py`
after the edit.

### 3.2 The O–C definition sentence — LEGITIMATE

**Proposed action.** Split Section 1.1 after "timing residuals,
conventionally called O–C". The definition of the observable should be its
own sentence.

### 3.3 "There is no corresponding way to make the star fire late" — LEGITIMATE

**Check.** Section 2.1 states it as a fact and, four sentences later, calls
the one-sided model an assumption rather than a result.

**Proposed action.** "In this picture there is no corresponding way to make
the star fire late."

### 3.4 "Thirty times the tunable beam" — LEGITIMATE

**Check.** The phrase is not defined. The Glashow absorbed fractions
(0.023, \(3\times 10^{-3}\), 0.064) are about thirty times the 1 EeV
deep-inelastic fractions on the same limbs (\(8\times 10^{-4}\),
\(10^{-4}\), \(2\times 10^{-3}\)).

**Proposed action.** "thirty times what the same limbs absorb from the
deep-inelastic beam at 1 EeV".

### 3.5 The column ratio stated, then the columns listed — LEGITIMATE

**Check.** The formula is the ratio
\(\Sigma_{\rm disk}/\Sigma_{\rm graze}\). The numbers that follow are the
grazing columns against the disk columns. The ratios are about
\(3\times 10^{3}\), \(3\times 10^{5}\) and \(10^{5}\). The disk columns are
not in Table 1, so dropping them would lose them.

**Proposed action.** State the three ratios, and keep the six columns in
parentheses after them. Leave the formula in the sentence; it is what the
ratios evaluate.

### 3.6 The Cepheid drift of \(+0.25\) s yr\(^{-1}\) — LEGITIMATE

**Check.** The sentence in Section 2.3 does not say the duty cycle. At
\(d = 1\) the same inputs give about \(0.5\) s yr\(^{-1}\). The \(0.25\)
is the \(d = 0.5\) case, which is the duty cycle used everywhere else in
Section 2.

**Proposed action.** "at a duty cycle of one half".

### 3.7 "How we weight what follows" — LEGITIMATE

**Check.** The last sentence of Section 2.7 says the mismatch shapes how we
weight what follows. The thing it should point at is that the null result
is carried by the RR Lyrae arm: that class is the one the geometry helps
least, and the only one the data can hear.

**Proposed action.** Replace the clause with "and it is why the null result
that follows is an RR Lyrae result".

### 3.8 The first pass, introduced mid-pipeline — LEGITIMATE

**Check.** "The first pass" appears without a definition in Section 4.3,
then again in Sections 4.9, 4.12, 6.1 and 6.2. Every count in the paper is
from the second run except where the text says otherwise (the two per cent
floor, the first-pass top ten, the catalog-flip control).

**Proposed action.** One sentence in the Section 4 preamble: the pipeline
was run twice; the earlier run is the first pass, and every count in this
paper is from the second unless the text says so.

### 3.9 "Which is why it drives at all" — LEGITIMATE, twice rather than three

**Check.** The phrase is in Section 2.2 ("that is why the zone drives at
all") and again in Section 2.4 ("which is why it drives at all"). A third
"drive" in the deposition paragraph ("heat leaks out within a cycle and
does not drive") is a different claim and should stay.

**Proposed action.** Keep the Section 2.2 sentence. In Section 2.4, replace
the clause with "(Section 2.2)".

### 3.10 The beat-control sentence, said twice — LEGITIMATE

**Check.** Section 4.9 ends the pipeline description with "none of 228 beat
injections into quiet TESS RR Lyrae stars is promoted by a decode route".
Section 4.12 states the same control in full, including that result.

**Proposed action.** Drop the parenthetical from Section 4.9 and leave a
pointer to Section 4.12. Keep the full statement there.

### 3.11 "Double or triple" the per-cycle scatter — LEGITIMATE, low corner is short of double

**Check.** A kick of 118–197 s adds white timing noise of \(\delta t/2\),
so 59–99 s, in quadrature with a 20–40 s scatter. The four corners:

| kick | scatter | factor |
|---:|---:|---:|
| 118 s | 40 s | 1.8 |
| 118 s | 20 s | 3.1 |
| 197 s | 40 s | 2.7 |
| 197 s | 20 s | 5.0 |

"Double or triple" is the middle of that range. "At least double" misses
the 1.8 corner.

**Proposed action.** "by a factor of about two to five".

### 3.12 TIC 147576267, scatter below the formal error — LEGITIMATE

**Check.** 375 s is 0.41 times 914 s. Section 4.4 inflates formal errors by
\(\sqrt{\chi^2_\nu}\) of the flux fit. The score parts for this star record
the same ratio ("O–C rms is 0.4× the timing error"). It is the only star
the paper quotes both numbers for. TIC 347799445, once Section 6.1 names
it, is a second case (108 s against 851 s), so the clause belongs on the
method, not as a special pleading for one star.

**Proposed action.** In the TIC 147576267 paragraph: the formal error is
the fit error inflated by \(\sqrt{\chi^2_\nu}\) (Section 4.4), and here
that factor sits above the scatter of the shifts.

---

## 4. Formatting and house style

### 4.1 Spelling — LEGITIMATE

**Check.** AAS house style is American. British forms still in the prose:
favoured, behaviour (several), colour and colour-coded, labelled, centre
(twice, "the star's centre"), centred, catalogued, modelling, analysed.
"Catalog" is already the dominant form. The keyword separator `---` and a
table's `Note---` stay; they are journal markup.

**Proposed action.** One pass to American spelling on those words. Do not
touch the markup.

### 4.2 Vetting class names — LEGITIMATE

**Check.** Classes are plain capitals everywhere except Section 4.9, where
KIC 5520878 is `\textsc{astro}`.

**Proposed action.** Change that one to ASTRO.

### 4.3 Table styles — LEGITIMATE, do not convert the environments

**Check.** Tables 1 and 2 are `tabular` because the preamble patches
AASTeX's `\tabular`, which shifts a centered table. Tables 3–6 are
`deluxetable`. The notes differ with the environments: a hand-set
`Note---` minipage against `\tablecomments`.

**Proposed action.** Leave Tables 1 and 2 as `tabular`. Set their note
minipages to the same size and width the deluxetable notes use, so the
notes match even though the rule weights will not.

### 4.4 Table 1, the ceiling column — LEGITIMATE

**Check.** The ceiling column equals the \(f_{\rm abs}\) column in every
row. The note already says the ceiling is \(f_{\rm abs}\).

**Proposed action.** Drop the column. The note keeps the sentence.

### 4.5 Table 5, the explanation cells — LEGITIMATE

**Proposed action.** "alias of the prewhitened 0.61 mode" and "subharmonic
of that alias, mode prewhitened", in place of "0.61 alias, mode
prewhitened" and "0.61 subharmonic, mode prewhitened".

### 4.6 Figure 1, \(H/R\approx 0.03\) — LEGITIMATE

**Check.** 0.03 is the RRab value. The Cepheid in the same section is
0.042, and the panel is generic.

**Proposed action.** "\(H/R\) of a few per cent would place it at the edge
of the disk".

### 4.7 Tarone, and Eddington and Plakidis — LEGITIMATE

**Check.** "Tarone's reduction" in Section 4.6 has no citation. The paper
is Tarone (1990), Biometrics 46, 515. It is not in either bib file.
`Eddington1929` already lists both authors. The first mention is
`\citet{Eddington1929}`; the later text says "the Eddington–Plakidis test".
If the bst renders the citation as both names, the surprise is already
gone. If it renders "Eddington (1929)", the sentence has to name Plakidis.

**Proposed action.** Add the Tarone entry and cite it at the reduction.
When that edit is made, validate it with `util/ads_validate.py` and
recompile. On the same compile, read the rendered Eddington citation: leave
it if both names appear, and write "Eddington and Plakidis (1929)" in the
prose if they do not.

### 4.8 Acknowledgments in the plural — LEGITIMATE

**Check.** The paper has one author. The paragraph says "the authors",
"hand-checked by the authors" and "the listed authors". It is this paper's
text, not a journal boilerplate that has to stay plural.

**Proposed action.** Make those three singular.

### 4.9 The named-word paragraph — LEGITIMATE

**Proposed action.** Split Section 4.6 at "The pass is also a statement
about named-word claims", so the across-star correction and the screening
decision are separate paragraphs.

### 4.10 The abstract's \(\beta\) — ALREADY THE RULE

**Check.** The abstract writes `$\beta$ Cephei` and uses no custom macros.
`util/check_abstract.py` reports no forbidden macros. That is required:
arXiv's metadata abstract does not expand `\newcommand`.

**Proposed action.** A one-line source comment above `\begin{abstract}`
saying the abstract is primitive TeX on purpose. No change to the text
beyond item 3.1.

---

## 5. Small things

### 5.1 "Remained largely a radio and, later, optical" — LEGITIMATE

**Proposed action.** "For its first half-century it was largely a radio and
then an optical enterprise."

### 5.2 "The catalog is all-sky; the search is not." — LEGITIMATE

**Check.** The line is already the last sentence of the paragraph. The
sentence before it ends on "absolute number", which is what blunts it.

**Proposed action.** Reword the Cepheid sentence so it does not end on
"number" ("The Cepheids are a larger fraction of their class and a much
smaller count."). Leave the punchline last.

### 5.3 "(6 for \(\beta\) Cephei)" — LEGITIMATE

**Proposed action.** "up to 12 harmonics, or 6 for \(\beta\) Cephei stars".

### 5.4 The sparse-timing sentence — LEGITIMATE

**Check.** The sentence that begins "First, each epoch's flux residual"
runs through the linearization, the bias on the rising branch, the block
average and the analysis, in one period.

**Proposed action.** End the main sentence at the block average. Give the
linearization bias its own sentence.

### 5.5 "Claim" for the jitter signature — LEGITIMATE

**Proposed action.** In the Section 6 opener, "signature" in place of
"claim" (both occurrences in that pair of sentences).

### 5.6 The QUIET score of 0.17 — LEGITIMATE, one more term

**Check.** TIC 147576267 has score 0.168. The parts are a quiet base of
0.06, a named-word bonus of 0.097, and \(0.30\times 0.287\times 0.137 =
0.012\) from the evidence term. Base plus bonus alone is 0.157, short of
the printed score.

**Proposed action.** "(quiet base 0.06, a named-word bonus of 0.10, and a
small evidence term)".

### 5.7 Conclusion 3 — LEGITIMATE

**Proposed action.** Split the sentence before "On the ground-based
surveys", so the injection result and the ground-based result are separate
sentences.

### 5.8 The cycle count in the closing paragraph — LEGITIMATE

**Check.** \(r\approx(3\times 30/0.7)^{2}\approx 1.5\times 10^{4}\) is the
only arithmetic done inside a prose sentence. The result (about
\(1.5\times 10^{4}\) repeats, \(4\times 10^{5}\) cycles, six centuries) is
what the paragraph needs.

**Proposed action.** Move the arithmetic to a footnote. Leave the result
in the sentence.
