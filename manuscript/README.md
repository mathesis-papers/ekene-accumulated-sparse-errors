# Accumulated sparse errors

Haobo Ma — ChronoAI Pte Ltd, Singapore  
Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore

## Results

The note addresses the two conjectures and the intervening pointwise claim in
Ekene E., *A bound on accumulated PSD contractions*, MathOverflow question
509794, revision 2, paragraphs 1–2, 2 April 2026.

The variable-correlation family has attained minima over arbitrary real
coefficients. It gives a 64-atom counterexample, the optimal pointwise envelope
`1/(1+sqrt(d))^2`, and accumulated finite-dimensional growth with sharp leading
constant one. A countably supported measure has a divergent accumulated series.
The classical pointwise estimate is credited to Maurey's empirical method,
Pisier (1981), Barron (1993), and the explicit greedy bounds of Barron, Cohen,
Dahmen and DeVore (2008). Exact references and source locators are in `note.tex`.

## Contents and build

- `note.tex`, `note.pdf`: statements and complete written proofs.
- `verify.py`, `verification.json`: exact rational checks and expected output.
- `metadata.json`: scientific metadata.
- `LICENSE`: CC BY 4.0 notice.
- `SHA256SUMS.json`: hashes and sizes of these files.

With Tectonic 0.17.0 installed:

```sh
tectonic note.tex
```

Alternatively, use a standard LaTeX distribution with the packages named in
`note.tex`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is included in the main source; no BibTeX step is needed.

## Exact arithmetic

```sh
python3 verify.py > regenerated-verification.json
cmp regenerated-verification.json verification.json
```

The half-correlation checks use dimensions 1–24, 32, 64 and 128 and every
sparsity level from 1 to the dimension. They check the harmonic total,
optimizer values and stationarity, signed square-completion identities and
full-atom representations. Nine normal-equation cases use rational Gaussian
elimination, including all 64 atoms.

The variable-correlation checks use dimensions 1–16, 32 and 64, every sparsity
level, and correlations `1/17`, `1/5`, `1/2`, `4/5`, `16/17`. They check the
shifted harmonic total, optimizers, signed square completion, and direct Gram
identities; 180 small systems are solved separately by rational elimination.
Additional shifted harmonic cases reach dimension 65536. Thirty rational
scalar-expectation controls include zero target energy and zero variance.

These finite checks support the displayed identities. The full parameter
range, all-real optimization, support correspondences, divergent countable
example and sharp asymptotic conclusions have written proofs in the note.

## Lean statements

The PDF links two fixed-version `.lean` files at commit
`88be2a120b67a9d9e724eaae970c8263d0a0e9a4`.
The quadratic-form declarations cover the half-correlation dimension-64
identity, attained all-real sparse minima, total exceeding one and the
negation of the finite assertion. The Hilbert-space declarations cover the
actual support, Bochner barycenter, equality of optimizing sets, exact minima,
zero tail, ordinary series total and negation of the Hilbert assertion.
The variable-correlation family and the other extensions have written proofs.

## Attribution and license

The posing source is licensed CC BY-SA 4.0; its statements are reformulated
and proved here without reproducing its prose. The note and accompanying
original verification code are licensed CC BY 4.0. See `LICENSE`.
