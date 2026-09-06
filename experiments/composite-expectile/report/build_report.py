#!/usr/bin/env python3
"""Build the standalone HTML report, embedding the figures as data URIs.

Run from the experiment directory:  python3 report/build_report.py
Writes report/report.html.
"""
import base64
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "out"


def fig(name, caption, alt):
    data = base64.b64encode((OUT / name).read_bytes()).decode()
    return (
        f'<figure class="wide">\n'
        f'  <img src="data:image/png;base64,{data}" alt="{alt}">\n'
        f'  <figcaption>{caption}</figcaption>\n'
        f'</figure>\n'
    )


HEAD = """<title>The Composite Expectile Estimator</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&display=swap">
<style>
:root {
  --ground:      #f4f7f8;
  --surface:     #ffffff;
  --surface-2:   #e9eff2;
  --ink:         #13202a;
  --ink-2:       #4b6373;
  --ink-3:       #7a8f9c;
  --rule:        #d5dee3;
  --rule-strong: #b3c2ca;
  --l1:          #1f5f8b;
  --l2:          #197a45;
  --ref:         #a53a2b;
  --flag:        #8a6d1f;
  --sans: "IBM Plex Sans", ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
  --serif: "Source Serif 4", Georgia, "Times New Roman", serif;
  --mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
  --measure: 68ch;
  --wide: 1120px;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --ground:      #0d141a;
    --surface:     #141d25;
    --surface-2:   #1a252e;
    --ink:         #e3ebf0;
    --ink-2:       #93a7b4;
    --ink-3:       #6d8290;
    --rule:        #253340;
    --rule-strong: #3a4c5a;
    --l1:          #6cb2e2;
    --l2:          #55c78a;
    --ref:         #e28d7e;
    --flag:        #d6b45c;
  }
}
:root[data-theme="dark"] {
  --ground:      #0d141a;
  --surface:     #141d25;
  --surface-2:   #1a252e;
  --ink:         #e3ebf0;
  --ink-2:       #93a7b4;
  --ink-3:       #6d8290;
  --rule:        #253340;
  --rule-strong: #3a4c5a;
  --l1:          #6cb2e2;
  --l2:          #55c78a;
  --ref:         #e28d7e;
  --flag:        #d6b45c;
}

* { box-sizing: border-box; }
body {
  background: var(--ground);
  color: var(--ink);
  font-family: var(--serif);
  font-size: 17px;
  line-height: 1.62;
  -webkit-font-smoothing: antialiased;
}
.page { max-width: var(--wide); margin: 0 auto; padding: 0 24px 96px; }

/* --- masthead ---------------------------------------------------------- */
header.masthead {
  border-bottom: 2px solid var(--ink);
  padding: 56px 0 22px;
  margin-bottom: 40px;
}
.eyebrow {
  font-family: var(--sans);
  font-size: 11.5px;
  font-weight: 600;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  color: var(--ink-3);
  margin-bottom: 14px;
}
h1 {
  font-family: var(--sans);
  font-weight: 600;
  font-size: clamp(32px, 5.2vw, 50px);
  line-height: 1.06;
  letter-spacing: -0.022em;
  margin: 0 0 16px;
  text-wrap: balance;
}
.standfirst {
  max-width: var(--measure);
  font-size: 19.5px;
  line-height: 1.55;
  color: var(--ink-2);
  margin: 0;
}
.byline {
  font-family: var(--sans);
  font-size: 13px;
  color: var(--ink-3);
  margin-top: 22px;
  display: flex;
  flex-wrap: wrap;
  gap: 6px 20px;
}

/* --- structure --------------------------------------------------------- */
main { display: flex; flex-direction: column; gap: 8px; }
section { scroll-margin-top: 20px; }
.col { max-width: var(--measure); }
h2 {
  font-family: var(--sans);
  font-weight: 600;
  font-size: 25px;
  letter-spacing: -0.012em;
  line-height: 1.2;
  margin: 56px 0 6px;
  text-wrap: balance;
  display: flex;
  gap: 14px;
  align-items: baseline;
}
h2 .num {
  font-family: var(--mono);
  font-size: 14px;
  font-weight: 500;
  color: var(--l1);
  flex: none;
  padding-top: 2px;
}
h3 {
  font-family: var(--sans);
  font-weight: 600;
  font-size: 17.5px;
  letter-spacing: -0.005em;
  margin: 34px 0 4px;
  color: var(--ink);
}
p { margin: 0 0 15px; }
h2 + p, h3 + p { margin-top: 10px; }
a { color: var(--l1); text-decoration-thickness: 1px; text-underline-offset: 2px; }
strong { font-weight: 600; }
code, .mono { font-family: var(--mono); font-size: 0.87em; }
p code, li code, td code {
  background: var(--surface-2);
  padding: 1px 4px;
  border-radius: 2px;
}
ul, ol { max-width: var(--measure); padding-left: 22px; margin: 0 0 15px; }
li { margin-bottom: 7px; }

/* --- display equations -------------------------------------------------- */
.eq {
  font-family: var(--mono);
  font-size: 14px;
  line-height: 1.75;
  background: var(--surface);
  border-left: 3px solid var(--rule-strong);
  padding: 16px 20px;
  margin: 20px 0;
  max-width: var(--measure);
  overflow-x: auto;
  white-space: pre;
  color: var(--ink);
}
.eq .hl { color: var(--ref); font-weight: 500; }

/* --- callouts ----------------------------------------------------------- */
.note {
  max-width: var(--measure);
  background: var(--surface);
  border: 1px solid var(--rule);
  border-left: 3px solid var(--flag);
  padding: 16px 20px;
  margin: 22px 0;
  font-size: 16px;
}
.note .lab {
  font-family: var(--sans);
  font-size: 11px;
  font-weight: 600;
  letter-spacing: 0.11em;
  text-transform: uppercase;
  color: var(--flag);
  display: block;
  margin-bottom: 6px;
}
.note p:last-child { margin-bottom: 0; }

/* --- verdict panel ------------------------------------------------------ */
.verdict {
  background: var(--surface);
  border: 1px solid var(--rule);
  border-top: 3px solid var(--l2);
  padding: 26px 28px;
  margin: 8px 0 12px;
}
.verdict h2 { margin-top: 0; font-size: 20px; }
.findings { display: grid; gap: 0; margin-top: 4px; }
.finding {
  display: grid;
  grid-template-columns: 88px 1fr;
  gap: 20px;
  padding: 15px 0;
  border-top: 1px solid var(--rule);
  align-items: start;
}
.finding:first-child { border-top: none; }
.finding .verdict-tag {
  font-family: var(--sans);
  font-size: 11px;
  font-weight: 600;
  letter-spacing: 0.07em;
  text-transform: uppercase;
  padding-top: 3px;
}
.t-yes { color: var(--l2); }
.t-no  { color: var(--ref); }
.t-mix { color: var(--flag); }
.finding p { margin: 0; font-size: 16.2px; }
@media (max-width: 620px) {
  .finding { grid-template-columns: 1fr; gap: 4px; }
}

/* --- tables ------------------------------------------------------------- */
.tablewrap { overflow-x: auto; margin: 22px 0; max-width: var(--wide); }
table {
  border-collapse: collapse;
  font-family: var(--sans);
  font-size: 13.5px;
  font-variant-numeric: tabular-nums;
  min-width: 100%;
}
caption {
  text-align: left;
  font-family: var(--sans);
  font-size: 13px;
  color: var(--ink-2);
  padding-bottom: 9px;
  max-width: var(--measure);
  line-height: 1.5;
}
th, td { padding: 7px 13px; text-align: right; white-space: nowrap; }
th:first-child, td:first-child { text-align: left; padding-left: 0; }
thead th {
  font-weight: 600;
  font-size: 12px;
  letter-spacing: 0.03em;
  color: var(--ink-2);
  border-bottom: 1.5px solid var(--rule-strong);
}
tbody td { border-bottom: 1px solid var(--rule); }
tbody tr:last-child td { border-bottom: 1.5px solid var(--rule-strong); }
tr.ref td { color: var(--ink-2); font-style: italic; }
td.win { color: var(--l2); font-weight: 600; }
td.lose { color: var(--ref); }
td.tie { font-weight: 600; }
.se { color: var(--ink-3); font-weight: 400; font-size: 0.88em; }

/* --- figures ------------------------------------------------------------ */
figure { margin: 30px 0; }
figure.wide { max-width: var(--wide); }
figure img {
  width: 100%;
  display: block;
  background: #ffffff;
  border: 1px solid var(--rule);
  border-radius: 2px;
}
figcaption {
  font-family: var(--sans);
  font-size: 13px;
  line-height: 1.55;
  color: var(--ink-2);
  margin-top: 10px;
  max-width: var(--measure);
}
figcaption b { color: var(--ink); font-weight: 600; }

/* --- footer ------------------------------------------------------------- */
footer {
  margin-top: 72px;
  padding-top: 22px;
  border-top: 1px solid var(--rule);
  font-family: var(--sans);
  font-size: 13px;
  color: var(--ink-3);
  max-width: var(--measure);
}
:focus-visible { outline: 2px solid var(--l1); outline-offset: 2px; }
@media (prefers-reduced-motion: reduce) { * { animation: none !important; transition: none !important; } }
</style>
"""

MASTHEAD = """
<div class="page">
<header class="masthead">
  <div class="eyebrow">Simulation study &middot; extreme value estimation</div>
  <h1>The Composite Expectile Estimator</h1>
  <p class="standfirst">Fitting a GEV by an expectile loss integrated against a tail-weighted
  measure, so that a misspecified body does not compromise the tail. Tested against maximum
  likelihood, L-moments, and the composite quantile estimator, on a Gaussian-contaminated GEV,
  over return periods from 2 to 1000 years.</p>
  <div class="byline">
    <span>2000 Monte Carlo replicates &middot; n = 50, 100, 200</span>
    <span>Truth: max(GEV, Normal)</span>
    <span>All code and outputs reproducible</span>
  </div>
</header>
<main>
"""

VERDICT = """
<section class="verdict">
<h2 style="display:block">What the experiment found</h2>
<div class="findings">
  <div class="finding"><div class="verdict-tag t-yes">Confirmed</div>
    <p><b>L2 fixes the variance explosion.</b> At matched weight, the composite
    <em>expectile</em> estimator's mean squared error for the 1000-year return level is
    <b>0.17&ndash;0.58&times;</b> that of the composite <em>quantile</em> estimator &mdash; a
    two- to six-fold improvement, at every weight tested. The L1 version's far-tail standard
    deviation is roughly double the L2 version's.</p></div>
  <div class="finding"><div class="verdict-tag t-yes">Confirmed</div>
    <p><b>Both remove the bias.</b> Maximum likelihood and L-moments underestimate the
    1000-year return level by 39% and 32% asymptotically under this contamination. Both
    composite estimators cut that to under 1%.</p></div>
  <div class="finding"><div class="verdict-tag t-mix">Qualified</div>
    <p><b>The weight was starting too late.</b> With the weight switching on near the median
    rather than at the 95th percentile, the composite expectile estimator <em>matches</em>
    the GEV MLE on MSE from T = 200 to T = 1000 (ratios 0.98 to 1.05, all within Monte Carlo
    error of 1) while carrying a quarter of its bias. In this contamination setting it never significantly beats it &mdash;
    but see below.</p></div>
  <div class="finding"><div class="verdict-tag t-no">New obstacle</div>
    <p><b>Expectiles are anchored to the mean</b>, so body contamination never fully leaves
    them. No weight function makes the L2 version asymptotically unbiased, whereas the L1
    version is exactly unbiased. And the raw L2 criterion is infinite once &xi; &ge; 1/2 &mdash;
    precisely the heavy-tailed regime the method exists to serve.</p></div>
  <div class="finding"><div class="verdict-tag t-yes">New</div>
    <p><b>Mixing the two losses helps.</b> A convex combination of the check functions elicits
    an intermediate functional &mdash; the &alpha;-elastile. Its contamination, bias and variance
    all interpolate smoothly, and an interior &alpha; beats <em>both</em> pure losses by 5 to 23%
    on MSE at every return period up to 500. The best &alpha; rises with the return period.</p></div>
  <div class="finding"><div class="verdict-tag t-yes">And</div>
    <p><b>Grafting onto an empirical body fixes the rest.</b> Used as the tail of a smooth graft
    with the empirical distribution as the body, the composite expectile fit goes from
    <b>1168&times;</b> the MLE's error at the 2-year level to <b>1.15&times;</b>, and beats the
    MLE across the whole tail (0.82&times; at T = 200, 0.84&times; at T = 1000, closer on 80% of
    datasets). Grafting the MLE onto the same body gains nothing, so it is the composite tail
    doing the work.</p></div>
  <div class="finding"><div class="verdict-tag t-yes">And</div>
    <p><b>There is a realistic setting where L2 wins outright.</b> When the body is
    near-degenerate &mdash; ordinary snowmelt maxima clustered tightly, a heavy rain-on-snow tail
    &mdash; the composite expectile estimator's typical error at the 1000-year level is
    <b>0.58&times;</b> the MLE's, and it is the closer of the two on <b>84%</b> of datasets. The
    composite quantile estimator is <em>worse</em> than the MLE there, winning 38%.</p></div>
</div>
</section>
"""

S1 = """
<section id="idea">
<h2><span class="num">01</span><span>The estimator, and a correction to the loss</span></h2>
<div class="col">
<p>The proposal is to fit a parametric family indexed by &theta; by minimising an expectile
loss integrated over levels, against a weight that is zero on the body of the distribution and
rises towards one in the tail:</p>
</div>
<div class="eq">S(&theta;) = &Sigma;&#8202;<sub>i</sub> &#8747;&#8202;<sub>0</sub><sup>1</sup>  w(p) <span class="hl">|</span>p &minus; I(y&#8202;<sub>i</sub> &lt; E(p|&theta;))<span class="hl">|</span> (y&#8202;<sub>i</sub> &minus; E(p|&theta;))<sup>2</sup>  dp</div>
<div class="col">
<p>where E(&middot;|&theta;) is the family's <em>expectile</em> function. The reasoning is that a
named parametric family is essentially always misspecified, and a global fitting criterion &mdash;
likelihood, L-moments &mdash; has to spread that misspecification across the whole distribution.
If only the tail matters, weighting the criterion towards the tail should buy a better tail fit.
The same construction with the pinball loss and the quantile function is the composite quantile
estimator; the question here is whether replacing L1 by L2 tames the variance explosion that made
the L1 version not worth using.</p>

<p>The absolute value marked in red is necessary and was missing from the proposal as stated.
Without it the factor (p &minus; I) is negative for every observation below E, so pushing E upward
sends the criterion to &minus;&infin;; minimising it numerically simply runs to whatever boundary
the optimiser is given. The ordering (y &minus; E) is correct as written &mdash; it is squared, so it
does not matter &mdash; and the indicator I(y &lt; E) is correct too. The pinball loss
(p &minus; I(y &lt; q))(y &minus; q) takes no absolute value and is already non-negative, which is
presumably where the slip came from.</p>

<p>With the absolute value, each integrand is the canonical strictly consistent scoring function
for its functional, so the integrated loss is minimised in expectation by a &theta; whose expectile
(respectively quantile) function matches the truth wherever w &gt; 0.</p>
</div>
</section>
"""

S2 = """
<section id="design">
<h2><span class="num">02</span><span>The simulation design</span></h2>
<div class="col">
<p>The truth is a Gaussian-contaminated GEV: the distribution of max(X, Z) with X ~ GEV and
Z ~ Normal independent, so its cdf is the product</p>
</div>
<div class="eq">F(x) = F&#8202;<sub>GEV</sub>(x; 0, 1, 0.2) &middot; &Phi;((x &minus; 1.5) / 0.8)</div>
<div class="col">
<p>The normal sits in the middle: it dominates the body and vanishes in the tail, so the GEV is
the <em>exactly correct tail model</em> and a <em>wrong body model</em>. That is the situation the
tail-focused idea is built for.</p>

<p>It is worth looking at what this actually does, because on a return-period axis the
contamination looks smaller than it is. On the quantile scale the two curves have converged by
about T = 15; the damage is all in the body, below the 90th percentile, where the truth's median
is 4.6&times; its GEV component's and its mean is 2.4&times;. That body is exactly what likelihood
and L-moments spend their fit on.</p>
</div>
FIG_DGP
<div class="col">
<p>The second row is the alternative setting suggested as a check, GEV(23, 22.2, 0.39) against
N(61.9, 16.9). It has a much heavier tail but, on the quantile scale, <em>milder</em>
contamination: 2.5% at T = 10 against 11% here, and a mean ratio of 1.6 against 2.4. The
setting used throughout is the more contaminated of the two.</p>

<h3>The weight function</h3>
<p>A smoothstep, zero below p&#8320; and one above p&#8321;. The primary design used
p&#8320; = 0.95 and p&#8321; = 0.98, placed where the <em>quantile</em> contamination dies out
(it leaves the quantiles 0.6% too high at p&#8320; and 0.00003% at p&#8321;). Section 6 shows
this choice, principled for L1, is badly wrong for L2, and section 5 sweeps it.</p>

<h3>Estimators and reporting</h3>
<p>GEV by maximum likelihood; GEV by L-moments (Hosking); composite quantile (L1); composite
expectile (L2). 2000 replicates at n = 50, 100 and 200, all four estimators fitted to identical
datasets, so every comparison below can be made <em>paired</em> &mdash; which matters, because the
sampling distribution of a far-tail return level is heavy-tailed and unpaired Monte Carlo error
would swamp the differences at issue. Every ratio quoted carries a paired Monte Carlo standard
error. The shape parameter is bounded to [&minus;0.45, 0.90] for all four estimators alike.</p>
</div>
</section>
"""

S3 = """
<section id="anchoring">
<h2><span class="num">03</span><span>Expectiles are anchored to the mean</span></h2>
<div class="col">
<p>This is the finding that most constrains the idea, and it is specific to L2.</p>

<p>The p-quantile of the contaminated truth equals the p-quantile of its GEV component as soon
as the normal's cdf reaches one. Quantiles are <em>tail-local</em>. Expectiles are not. The
p-expectile solves</p>
</div>
<div class="eq">k &middot; &phi;(x) + m &minus; x = 0,   k = (2p &minus; 1)/(1 &minus; p),   &phi;(x) = E[(Y &minus; x)<sup>+</sup>]</div>
<div class="col">
<p>and while &phi;(x) for large x sees only the tail, the mean m sees everything &mdash; including
the contaminated body. So the truth's expectiles differ from its GEV component's at <em>every</em>
level, by an amount that decays only like (1 &minus; p)<sup>&xi;</sup>.</p>
</div>
<div class="tablewrap">
<table>
<caption>How much of the body contamination is still present in each functional. The quantile
discrepancy reaches machine precision by the 50-year level; the expectile discrepancy is still
9.7% there, and 1.4% at 10&#8202;000 years.</caption>
<thead><tr><th>level p</th><th>return period</th><th>quantile still off by</th><th>expectile still off by</th></tr></thead>
<tbody>
<tr><td>0.95</td><td>20 yr</td><td>0.56%</td><td class="lose">16.0%</td></tr>
<tr><td>0.98</td><td>50 yr</td><td>0.00003%</td><td class="lose">9.7%</td></tr>
<tr><td>0.99</td><td>100 yr</td><td>&asymp; 0</td><td class="lose">7.0%</td></tr>
<tr><td>0.999</td><td>1000 yr</td><td>&asymp; 0</td><td class="lose">2.9%</td></tr>
<tr><td>0.9999</td><td>10&#8202;000 yr</td><td>&asymp; 0</td><td class="lose">1.4%</td></tr>
</tbody></table>
</div>
<div class="col">
<p>The consequence is structural: <strong>there is no weight function that makes the composite
expectile estimator asymptotically unbiased under body contamination</strong>, whereas the
composite quantile estimator is exactly unbiased. Fitting each estimator to a single sample of
2&times;10<sup>6</sup> gives the parameter each converges to, and hence the bias no sample size
removes:</p>
</div>
<div class="tablewrap">
<table>
<caption>Asymptotic targets and the return-level bias they imply, at the primary weight
(p&#8320; = 0.95). The tail-focused idea works exactly as intended &mdash; both composite
estimators essentially eliminate the 20&ndash;40% bias that the global methods carry.</caption>
<thead><tr><th>estimator</th><th>&mu;</th><th>&sigma;</th><th>&xi;</th><th>bias at T = 100</th><th>at T = 1000</th></tr></thead>
<tbody>
<tr><td>GEV maximum likelihood</td><td>1.397</td><td>0.890</td><td>0.062</td><td class="lose">&minus;18.7%</td><td class="lose">&minus;39.1%</td></tr>
<tr><td>GEV L-moments</td><td>1.384</td><td>0.819</td><td>0.120</td><td class="lose">&minus;15.1%</td><td class="lose">&minus;31.6%</td></tr>
<tr><td>Composite quantile (L1)</td><td>0.011</td><td>0.999</td><td>0.200</td><td class="win">+0.02%</td><td class="win">&minus;0.08%</td></tr>
<tr><td>Composite expectile (L2)</td><td>0.832</td><td>0.879</td><td>0.219</td><td>+3.4%</td><td>+0.7%</td></tr>
<tr class="ref"><td>true GEV tail</td><td>0</td><td>1</td><td>0.20</td><td>&mdash;</td><td>&mdash;</td></tr>
</tbody></table>
</div>
<div class="col">
<p>L1 recovers the true tail parameters to three decimals. L2 lands a few percent away, and that
gap is the mean-anchoring above.</p>
</div>
</section>
"""

S4 = """
<section id="headline">
<h2><span class="num">04</span><span>Return levels, and the price of unbiasedness</span></h2>
<div class="col">
<p>Bias is only half the account. The panels below show the sampling distribution of each
estimator's return-level curve at n = 100, under the primary weight.</p>
</div>
FIG_RL
<div class="col">
<p>Maximum likelihood has a tight band whose centre flattens badly below the truth &mdash; the
20&ndash;40% bias of the previous section, visible as a curve that stops rising. The composite
estimators track the truth far better in the median and pay for it in width. The two lower-right
panels use a weight capped so that it does not extend past the largest observation (section 6);
the capped L1 median follows the truth almost exactly to T = 1000, with a 95% band running from
5 to beyond 30.</p>
</div>
FIG_MSE
<div class="col">
<p>On mean squared error, at this weight, no composite design beats maximum likelihood anywhere
from T = 20 to T = 1000. The trade each makes against the MLE at T = 1000, in MSE units, is the
compact way to see why:</p>
</div>
<div class="tablewrap">
<table>
<caption>Bias&ndash;variance accounting for the 1000-year return level, n = 100.
The true return level is 14.90.</caption>
<thead><tr><th>estimator</th><th>bias&sup2;</th><th>variance</th><th>MSE</th><th>bias&sup2; bought</th><th>variance paid</th></tr></thead>
<tbody>
<tr><td>GEV maximum likelihood</td><td>33.5</td><td>8.2</td><td>41.7</td><td>&mdash;</td><td>&mdash;</td></tr>
<tr><td>Composite quantile (L1)</td><td>3.2</td><td>76.1</td><td>79.3</td><td class="win">&minus;30.3</td><td class="lose">+67.9</td></tr>
<tr><td>Composite expectile (L2)</td><td>12.1</td><td>33.7</td><td>45.8</td><td class="win">&minus;21.4</td><td>+25.5</td></tr>
</tbody></table>
</div>
<div class="col">
<p>L1 pays 2.2 units of variance for every unit of squared bias it removes. L2 pays 1.2.
<strong>Switching to L2 turns a badly losing trade into a nearly even one</strong> &mdash; which is
the direction the L1-to-L2 argument predicts, and about half of the improvement needed.</p>
</div>
</section>
"""

S5 = """
<section id="tailfocus">
<h2><span class="num">05</span><span>How much is the tail focus, and how much is L1 versus L2?</span></h2>
<div class="col">
<p>The previous section fixes the weight at one setting, which conflates two things: the choice
of functional and the strength of the tail focus. Sweeping the lower edge p&#8320; from 0
(w &equiv; 1, no tail focus at all) up to 0.95 separates them. At p&#8320; = 0 the composite
quantile estimator is exactly CRPS minimisation, which is the like-for-like benchmark against
L-moments: both are global, quantile-flavoured, and neither has any tail focus.</p>
</div>
FIG_TF
<div class="col">
<p>The two panels behave completely differently in the far tail. Every L2 curve converges to
roughly 1.05 at T = 1000 whatever the weight; every L1 curve diverges upward, the worst of them
to 6&times; the MLE. That contrast <em>is</em> the answer to the original question.</p>
</div>
<div class="tablewrap">
<table>
<caption>Paired comparison of the two losses at matched weight, same 2000 datasets:
MSE(L2) / MSE(L1) for the 1000-year return level, with the Monte Carlo standard error of the
paired difference. L2 wins decisively at every weight.</caption>
<thead><tr><th>weight starts at</th><th>MSE(L2) / MSE(L1) at T = 1000</th><th>sd(L2)</th><th>sd(L1)</th></tr></thead>
<tbody>
<tr><td>p&#8320; = 0 &nbsp;(w &equiv; 1)</td><td class="win">0.83 <span class="se">&plusmn; 0.04</span></td><td>6.72</td><td>4.32</td></tr>
<tr><td>p&#8320; = 0.50</td><td class="win">0.23 <span class="se">&plusmn; 0.07</span></td><td>6.49</td><td>13.51</td></tr>
<tr><td>p&#8320; = 0.70</td><td class="win">0.17 <span class="se">&plusmn; 0.06</span></td><td>6.32</td><td>15.03</td></tr>
<tr><td>p&#8320; = 0.80</td><td class="win">0.21 <span class="se">&plusmn; 0.06</span></td><td>6.26</td><td>13.80</td></tr>
<tr><td>p&#8320; = 0.90</td><td class="win">0.38 <span class="se">&plusmn; 0.06</span></td><td>6.15</td><td>11.05</td></tr>
<tr><td>p&#8320; = 0.95</td><td class="win">0.58 <span class="se">&plusmn; 0.05</span></td><td>5.81</td><td>8.73</td></tr>
</tbody></table>
</div>
<div class="col">
<p>Two things stand out. First, L2's far-tail standard deviation is <em>almost invariant</em> to
the weight (5.8 to 6.7 across the whole sweep) while L1's swings between 4.3 and 15.0. The
expectile criterion is simply a much better conditioned thing to minimise under tail weighting,
because expectiles are smooth functionals of the whole sample while quantiles lean on a handful
of order statistics.</p>

<p>Second, the sign of the comparison flips at p&#8320; = 0. With no tail focus, the
<em>quantile</em> version has the lower variance (ratio 0.83 in L2's favour on MSE, but sd 4.32
against 6.72 &mdash; L2 wins on MSE there only through bias). So the L2 advantage is not a general
property of expectiles; it is specific to tail-weighted criteria, which is exactly where the
proposal lives.</p>

<h3>Where the weight should start</h3>
<p>The primary design put p&#8320; at 0.95, where the quantile contamination dies out. That is
correct reasoning for L1 and wrong for L2 &mdash; section 3 showed the expectile contamination is
still 16% there, so restricting to the extreme tail costs a great deal of data and buys almost no
bias reduction. The sweep bears that out.</p>
</div>
<div class="tablewrap">
<table>
<caption>MSE relative to GEV MLE for the 1000-year return level, n = 100, with the Monte Carlo
standard error of the paired difference. Ratios within about two standard errors of 1 are ties.
This is a post-hoc sweep, not a pre-registered choice &mdash; read it as a sensitivity analysis.</caption>
<thead><tr><th>weight starts at</th><th>composite quantile (L1)</th><th>composite expectile (L2)</th><th>median &xi;&#770; (L2)</th></tr></thead>
<tbody>
<tr><td>p&#8320; = 0 &nbsp;(w &equiv; 1)</td><td class="lose">1.43 <span class="se">&plusmn; 0.02</span></td><td class="lose">1.19 <span class="se">&plusmn; 0.05</span></td><td>0.162</td></tr>
<tr><td>p&#8320; = 0.50</td><td class="lose">4.64 <span class="se">&plusmn; 0.36</span></td><td class="tie">1.05 <span class="se">&plusmn; 0.05</span></td><td class="win">0.209</td></tr>
<tr><td>p&#8320; = 0.70</td><td class="lose">6.31 <span class="se">&plusmn; 0.39</span></td><td class="tie">1.04 <span class="se">&plusmn; 0.05</span></td><td>0.164</td></tr>
<tr><td>p&#8320; = 0.80</td><td class="lose">5.11 <span class="se">&plusmn; 0.33</span></td><td class="tie">1.06 <span class="se">&plusmn; 0.05</span></td><td>0.106</td></tr>
<tr><td>p&#8320; = 0.90</td><td class="lose">2.94 <span class="se">&plusmn; 0.22</span></td><td>1.11 <span class="se">&plusmn; 0.05</span></td><td>&minus;0.029</td></tr>
<tr><td>p&#8320; = 0.95</td><td class="lose">1.90 <span class="se">&plusmn; 0.13</span></td><td>1.10 <span class="se">&plusmn; 0.04</span></td><td>&minus;0.196</td></tr>
<tr class="ref"><td>GEV L-moments</td><td colspan="2">1.11</td><td>0.095</td></tr>
<tr class="ref"><td>GEV maximum likelihood</td><td colspan="2">1.00 (reference)</td><td>0.054</td></tr>
</tbody></table>
</div>
<div class="col">
<p>With the weight switching on near the median, the composite expectile estimator becomes
<strong>indistinguishable from the GEV MLE on MSE across the whole far tail</strong>. Following
the p&#8320; = 0.50 fit across return periods, the ratio to the MLE is 1.07 &plusmn; 0.04 at
T = 100, 0.98 &plusmn; 0.04 at T = 200, 0.98 &plusmn; 0.04 at T = 500 and 1.05 &plusmn; 0.05 at
T = 1000. It never significantly beats it. But it ties it while carrying a quarter of its bias
(&minus;1.39 against &minus;5.79 at T = 1000) and recovering the shape parameter almost exactly:
median &xi;&#770; = 0.209 against a true 0.20, where the MLE returns 0.054 and L-moments 0.095.</p>

<div class="note"><span class="lab">Worth noting</span>
<p>An estimator that ties on MSE but recovers &xi; correctly is not the same product as one that
wins on MSE. If what a study needs is a defensible shape parameter and a return-level curve that
does not flatten, that is a different and possibly more useful thing than a marginally smaller
squared error.</p></div>
</div>
</section>
"""

S6 = """
<section id="weights">
<h2><span class="num">07</span><span>Which weight functions are admissible</span></h2>
<div class="col">
<p>The weight cannot be chosen freely, and there are two separate conditions &mdash; one on the
weight, one on the distribution being fitted.</p>

<h3>The condition on the weight</h3>
<p>Writing u = 1 &minus; p, the population risk at level p decays as</p>
</div>
<div class="eq">quantile  (L1):   r&#8202;<sub>p</sub> ~ u<sup>1 &minus; &xi;</sup>
expectile (L2):   r&#8202;<sub>p</sub> ~ u<sup>1 &minus; 2&xi;</sup></div>
<div class="col">
<p>with the measured log-log slopes converging to these to three decimals as u &rarr; 0. So for a
weight that blows up like w(p) ~ (1 &minus; p)<sup>&minus;a</sup>, the criterion is finite if and
only if</p>
</div>
<div class="eq">quantile  (L1):   a &lt; 2 &minus; &xi;
expectile (L2):   a &lt; 2 &minus; 2&xi;</div>
<div class="col">
<p>At &xi; = 0.2 that is a &lt; 1.8 and a &lt; 1.6, so a weight going like 1/(1 &minus; p) is
admissible for both and 1/(1 &minus; p)&sup2; for neither. Note the L2 threshold is the tighter of
the two. Taking the worst case over the range where the quantile version is defined at all,
&xi; &uarr; 1, gives a &lt; 1 &mdash; which reproduces the L1 result reported in the original
thesis work.</p>

<p>The weight used throughout this study is bounded &mdash; identically 0 below p&#8320; and
identically 1 above p&#8321; &mdash; so it clears both conditions with room to spare at every
shape the estimators visit. Nothing here requires a weight bounded by 1; that is simply what was
used.</p>

<h3>The condition on the distribution &mdash; and this one separates L1 from L2</h3>
<p>The quantile risk needs only a finite first moment, &xi; &lt; 1. The expectile risk needs a
finite <em>second</em> moment, &xi; &lt; 1/2, because E[(e &minus; Y)&sup2; 1{Y &lt; e}] is
infinite for every e once the variance is:</p>
</div>
<div class="tablewrap">
<table>
<caption>Population risk at p = 0.99 for a GEV truth, by shape parameter.</caption>
<thead><tr><th>&xi;</th><th>quantile risk</th><th>expectile risk</th></tr></thead>
<tbody>
<tr><td>0.30</td><td>0.146</td><td>1.914</td></tr>
<tr><td>0.45</td><td>0.285</td><td>24.799</td></tr>
<tr><td>0.50</td><td>0.364</td><td class="lose">&infin;</td></tr>
<tr><td>0.60</td><td>0.623</td><td class="lose">&infin;</td></tr>
</tbody></table>
</div>
<div class="col">
<p>No weight function repairs this &mdash; the criterion is +&infin; at every &theta;. So the L2
version, taken literally, is undefined exactly in the heavy-tailed regime the method exists to
serve. For flood frequency that is not a remote corner of the parameter space: the alternative
setting examined here has &xi; = 0.39.</p>

<p>The repair is standard: minimise the loss <em>difference</em> against a fixed reference
function,</p>
</div>
<div class="eq">&Sigma;&#8202;<sub>i</sub> &#8747; w(p) [ &rho;&#8202;<sub>p</sub>(y&#8202;<sub>i</sub>, T(p|&theta;)) &minus; &rho;&#8202;<sub>p</sub>(y&#8202;<sub>i</sub>, T&#8202;<sub>ref</sub>(p)) ] dp</div>
<div class="col">
<p>which has the same minimiser and is finite whenever E|Y| &lt; &infin;. In this study &xi; = 0.2,
so the raw criterion is finite and nothing was needed &mdash; but it should be built in before the
estimator meets anything heavier.</p>

<h3>A third trap: weight beyond the data</h3>
<p>A weight left at one all the way to p = 1 puts a large share of its mass at levels
<em>above the largest observation</em>. Out there the empirical quantile and expectile functions
have both saturated at the sample maximum, so the loss can only pull the fitted tail down. The
fitted shape parameter tracks that share almost exactly:</p>
</div>
<div class="tablewrap">
<table>
<caption>Median fitted shape against the share of the weight sitting above the largest observation,
for the p&#8320; = 0.95 design. True &xi; = 0.2. Both columns converge to it, so this is
finite-sample behaviour, not a defect of the estimator.</caption>
<thead><tr><th>n</th><th>weight above the largest observation</th><th>median &xi;&#770;, L1</th><th>median &xi;&#770;, L2</th></tr></thead>
<tbody>
<tr><td>50</td><td class="lose">45%</td><td class="lose">&minus;0.82</td><td class="lose">&minus;0.52</td></tr>
<tr><td>100</td><td class="lose">27%</td><td>&minus;0.24</td><td>&minus;0.19</td></tr>
<tr><td>200</td><td>13%</td><td>&minus;0.01</td><td>&minus;0.06</td></tr>
<tr><td>500</td><td>5%</td><td>0.17</td><td>0.09</td></tr>
<tr><td>10&#8202;000</td><td>0.3%</td><td class="win">0.21</td><td class="win">0.21</td></tr>
</tbody></table>
</div>
<div class="col">
<p>At n = 50 this pins 89% of L1 fits and 56% of L2 fits onto the shape bound. Capping the weight
at p = 1 &minus; 1/n fixes the shape (median &xi;&#770; moves from &minus;0.26 to 0.19 for L1 and
&minus;0.20 to 0.06 for L2) &mdash; but for L1 it then <em>worsens</em> the far-tail MSE by an
order of magnitude, because the downward drag had been acting as accidental shrinkage. L2 survives
capping intact. Whatever else is done with the estimator, <strong>the weight should not extend
past where the data can speak</strong>.</p>
</div>
FIG_SHAPE
</section>
"""

S8 = """
<section id="recommend">
<h2><span class="num">10</span><span>What to do with this</span></h2>
<div class="col">
<p><strong>The obstacle is variance, not bias.</strong> Both composite estimators solve the bias
problem completely &mdash; against a body-contaminated truth they remove a 39% underestimate of
the 1000-year level that maximum likelihood cannot. What neither solves on its own is the variance
they pay for it.</p>

<p><strong>Use L2, not L1, whenever the criterion is tail-weighted.</strong> This is the clearest
result in the study: two- to six-fold lower far-tail MSE at every weight tested, and a standard
deviation that barely moves as the weight changes. In the near-degenerate-body setting it is the
difference between beating the MLE on 84% of datasets and losing on 62% of them.</p>

<p><strong>Put the weight much lower than the quantile reasoning suggests.</strong> Because
expectile contamination decays like (1 &minus; p)<sup>&xi;</sup> rather than vanishing, restricting
an expectile criterion to the extreme tail costs data without buying bias reduction. Here the
useful range was p&#8320; &asymp; 0.5&ndash;0.8, not 0.95. This is not a detail: it is worth more
than any other single choice in the study.</p>

<p><strong>Mix the losses, with &alpha; increasing in the return period.</strong> At a
well-placed weight the interior optimum is real and worth 5&ndash;23%. Roughly &alpha; &asymp; 0.2
for T near the record length, rising to pure L2 for T an order of magnitude beyond it.</p>

<p><strong>Cap the weight at the level of the largest observation</strong>, and
<strong>build in the loss-difference form before going near heavy tails</strong>, where the raw L2
criterion is infinite for &xi; &ge; 1/2.</p>

<h3>Where the method earns its keep</h3>
<p>The contest is governed by
(&xi;<sub>true</sub> &minus; &xi;<sub>MLE</sub>) / sd(&xi;&#770;<sub>composite</sub>), not by the
size of the global fit's bias alone. So the method pays off when the truth drags a global fit's
shape hard <em>while leaving a large clean region for the weighted criterion to read</em> &mdash;
a sharply concentrated body with a heavy tail, which is a recognisable hydrological situation.
It does <em>not</em> pay off simply because the tail is heavier: heavier tails inflate the
composite estimators' variance at the same rate they inflate the MLE's bias, which is why the
grafted setting of section 6 favoured the MLE despite a 55% bias.</p>

<h3>What looks most promising next</h3>
<ul>
<li><strong>Reduce the effective number of free parameters over the weighted region.</strong> The
loss surface is a near-flat ridge in (&mu;, &sigma;, &xi;): six multi-starts agree to five decimals
on a criterion that varies only in the fourth decimal as &xi; ranges from &minus;0.15 to +0.45.
Almost all of the variance is that ridge, and attacking it directly will do more than any further
loss-function tinkering.</li>
<li><strong>Correct the mean anchoring.</strong> The L2 bias is entirely the gap between the
truth's mean and the model's. Matching the model's expectile function against a mean-corrected
empirical target would remove L2's one structural disadvantage against L1.</li>
<li><strong>Report the shape parameter, not only the MSE.</strong> At p&#8320; = 0.5 the expectile
fit recovers &xi; = 0.209 against a true 0.20 where likelihood returns 0.054. On the question
"does the fitted tail have the right heaviness", the two are not close, and an estimator that ties
on squared error while getting &xi; right is a different and arguably more useful product.</li>
<li><strong>Choose the loss criterion deliberately.</strong> Squared error at a 1000-year return
level is dominated by a handful of replicates. If an occasional overestimate of a design flood is
less costly than a systematic underestimate of every one &mdash; which in flood frequency it
usually is &mdash; then median absolute error or a head-to-head win rate is the criterion to
optimise, and by those the case for the method is much stronger than by MSE.</li>
</ul>
</div>
</section>
"""

S9 = """
<section id="methods">
<h2><span class="num">11</span><span>Methods and reproducibility</span></h2>
<div class="col">
<h3>Machinery</h3>
<p>The expectile function of the GEV is needed on a grid of levels inside an optimiser, roughly
10<sup>8</sup> times across the study, so it is computed in closed form rather than numerically.
Substituting z = (1 + &xi;(t &minus; &mu;)/&sigma;)<sup>&minus;1/&xi;</sup> into the partial moment
and integrating by parts gives, for &xi; &lt; 1 and &xi; &ne; 0,</p>
</div>
<div class="eq">&phi;(x) = (&sigma;/&xi;) [ &Gamma;<sub>lower</sub>(1 &minus; &xi;, z<sub>x</sub>) &minus; (1 &minus; e<sup>&minus;z<sub>x</sub></sup>) z<sub>x</sub><sup>&minus;&xi;</sup> ]</div>
<div class="col">
<p>with &Gamma;<sub>lower</sub>(s, a) = pgamma(a, s)&middot;&Gamma;(s) &mdash; one line in R, exact
to machine precision, and about 90&times; faster than integrating the survival function. The same
technique gives the second-order partial moment used to decide when the losses are finite. The
expectile itself is then a safeguarded Newton iteration on the identification equation, ported to
C++ for the simulation; the R reference implementation is retained and checked against it.</p>

<p>The composite losses are evaluated by cumulative sums over the sorted sample rather than an
n&times;K matrix, which makes them O(K log n) and identical to the direct definition to
10<sup>&minus;15</sup>. Levels are integrated by a composite Gauss-Legendre rule &mdash; 24 panels
of 8 nodes &mdash; under the substitution 1 &minus; p = (1 &minus; p&#8320;)s&sup2;, which absorbs
the unbounded derivative at p = 1 and localises to one panel each the kinks where the fitted
functional crosses an observation. Refining the rule further moves a fitted 100-year return level
by under 0.15%.</p>

<h3>Checks</h3>
<p>Every component is validated in a re-runnable script: the closed-form partial moments against
independent numerical integration and against endpoint and derivative identities; the expectile
solver against its own identification equation, against the R/C++ pair, and against the
<code>distionary</code> package; the fast losses against their direct definitions; the quadrature
against refinement; the optimiser against multi-start; and the simulated truth against its
analytic quantiles and mean.</p>

<p>Optimiser reliability was a live concern, since the composite loss surface is nearly flat. It
is not the explanation for anything reported here: five widely separated starts converge to the
same point to five decimals for both losses.</p>

<h3>A bug in the reference implementation</h3>
<p>Cross-checking turned up a genuine fault in the <code>expectiles</code> branch of
<code>probaverse/distionary</code>. For &xi; &gtrsim; 0.5 it returns silently wrong expectiles
&mdash; for GEV(0, 1, 0.7) at p = 0.7 it gives 5.690 where the answer is 5.129, and the returned
values are exactly 2, 8 and 32 times the mean, which are bracket-expansion endpoints rather than
roots. The cause is two things combining: <code>integrate(survival, x, Inf)</code> fails
sporadically on a heavy tail and returns NaN for particular lower limits while neighbouring ones
are fine; and the bracketing routine advances on <code>if (g(hi) &lt;= 0.0) return true;</code>,
where <code>NaN &lt;= 0</code> is false, so a single failed integration makes the search step over
the root. Recommended fixes: treat a non-finite g as a failure rather than letting it drive the
search, and replace the numerical partial moment with the closed form above.</p>
</div>
</section>
"""

FOOTER = """
</main>
<footer>
<p>Simulation study for a research sabbatical postdoc project on rain-on-snow flooding. Nothing
here depends on that analysis &mdash; the experiment is self-contained. Code, outputs and the full
validation log live in <code>experiments/composite-expectile/</code>.</p>
</footer>
</div>
"""


def build():
    body = [
        HEAD, MASTHEAD, VERDICT, S1,
        S2.replace("FIG_DGP", fig(
            "fig-dgp.png",
            "<b>The data-generating process.</b> Left, densities of the two components and of "
            "their maximum; right, the frequency&ndash;magnitude curves of all three. Top row, "
            "the setting used throughout; bottom, the heavier alternative. In both, the truth "
            "and its GEV component are indistinguishable beyond about T = 15, while the "
            "densities differ completely &mdash; the contamination is real, and it all lives in "
            "the body.",
            "Densities and frequency-magnitude curves for two contaminated-GEV settings")),
        S3,
        S4.replace("FIG_RL", fig(
            "fig-headline-return-level.png",
            "<b>Return period against return level, n = 100.</b> Solid black is the truth, the "
            "coloured line the median estimate, the band the central 95% of 2000 replicates. "
            "The two global methods are tight and bent; the composite estimators are centred and "
            "wide.",
            "Return level curves with 95 percent bands for six estimators"))
         .replace("FIG_MSE", fig(
            "fig-headline-mse.png",
            "<b>Mean squared error over return period.</b> Left, absolute MSE; right, MSE "
            "relative to the GEV MLE, where below 1 would beat it. At this weight, nothing does "
            "beyond T = 20.",
            "MSE and MSE ratio against return period")),
        S5.replace("FIG_TF", fig(
            "fig-tail-focus.png",
            "<b>Sweeping the tail focus.</b> MSE relative to the GEV MLE as the weight's lower "
            "edge moves from p&#8320; = 0 (w &equiv; 1, palest) to 0.95 (darkest); dashed orange "
            "is GEV L-moments. Every expectile curve converges to about 1.05 in the far tail "
            "whatever the weight; every quantile curve diverges upward.",
            "MSE ratio against return period for six weight settings, quantile and expectile")),
        S6.replace("FIG_SHAPE", fig(
            "fig-shape.png",
            "<b>Fitted shape against sample size</b>, with a wider shape bound than the main run "
            "so the bound itself cannot be the explanation. Both losses converge to the true "
            "&xi; = 0.2; the small-sample bias is the extrapolation drag described above.",
            "Boxplots of the fitted shape parameter against sample size")),
        S6B.replace("FIG_SHARP", fig(
            "fig-sharp.png",
            "<b>A near-degenerate body.</b> Left, the truth: a tight normal at 2.5 with a heavy "
            "GEV tail. Middle, MSE relative to the GEV MLE for both losses at three weights. "
            "Right, the median fitted return-level curve &mdash; maximum likelihood flattens, the "
            "composite expectile fit tracks the truth.",
            "Density, MSE ratio and median return level for the sharp-contrast setting")),
        S7E.replace("FIG_EL", fig(
            "fig-elastile.png",
            "<b>The &alpha;-sweep.</b> Top row, the well-placed weight (p&#8320; = 0.5); bottom "
            "row, the original tail-focused one (p&#8320; = 0.95). Left, MSE against return "
            "period for each &alpha;; middle and right, MSE against &alpha; at T = 100 and "
            "T = 1000, with the GEV MLE and L-moments marked. The interior minimum at "
            "p&#8320; = 0.5 is the elastic-net effect; at p&#8320; = 0.95 it is largely absent.",
            "MSE against return period and against alpha for the elastile at two weights")),
        SG.replace("FIG_SG", fig(
            "fig-smooth-graft.png",
            "<b>Grafting repairs the body and improves the tail.</b> Dashed lines are ungrafted "
            "composite fits, solid lines the same fits used as the tail of a smooth graft with an "
            "empirical body. Left, MSE relative to the GEV MLE on a log scale &mdash; note the "
            "ungrafted p&#8320;=0.95 expectile fit starting above 1000 at T = 2 and the grafted "
            "version sitting at 1.15. Right, median absolute error. The red line is the control: "
            "the MLE grafted onto the same body, which gains nothing.",
            "MSE ratio and median error ratio for grafted and ungrafted composite fits")),
        S8, S9, FOOTER,
    ]
    return "".join(body)



S6B = """
<section id="whenwins">
<h2><span class="num">06</span><span>When is the family wrong enough for this to win?</span></h2>
<div class="col">
<p>Everything so far says the tail-focused estimators tie the MLE at best. The natural next
question is whether that is a property of the method or of the test: is there a realistic truth
whose tail bias is severe enough that a global fit <em>must</em> lose? Two attempts, one that
failed instructively and one that worked.</p>

<h3>Making the family more wrong does not help by itself</h3>
<p>The first attempt replaces max(GEV, Normal) with a <em>graft</em>: above a join level the truth
IS the GEV, exactly, and below it the body is a different distribution carrying no information
about the tail at all. That removes the leak in the max construction, where the GEV component
still shapes the upper body and so hands a global fit some tail information for free.
Hydrologically it is the standard mixed-population picture, and it fits the rain-on-snow story
directly: in most years the annual maximum is an ordinary snowmelt peak from a much less variable
process, and in the rest it is a rain-on-snow event from a genuinely heavy-tailed one.</p>

<p>With GEV(0, 1, 0.4) above the 90th percentile and a truncated normal below, maximum likelihood
is 55% low at the 1000-year level &mdash; far worse than the 39% in the main setting. And every
composite estimator does <em>worse</em> against it, not better: MSE ratios of 2.8 to 4.5 at
T = 1000. The reason is that the heavier tail which inflates the MLE's bias inflates the composite
estimators' variance at the same time. Their standard deviation at T = 1000 is 35 to 45 on a true
return level of 37.</p>

<div class="note"><span class="lab">The governing ratio</span>
<p>To first order the fractional return-level error at large T behaves like log(T) times the error
in &xi;. So the contest is decided not by the size of the MLE's bias but by
<b>(&xi;<sub>true</sub> &minus; &xi;<sub>MLE</sub>) / sd(&xi;&#770;<sub>composite</sub>)</b> &mdash;
and that ratio is roughly constant in T, which is exactly why the MSE-ratio curves flatten out in
the far tail. Heavier tails move both terms together. What is needed is a truth that drags the
global fit's shape hard <em>while leaving a large, clean region for the weighted criterion to
read</em>.</p></div>

<h3>A near-degenerate body does</h3>
<p>The second attempt keeps the max construction and simply makes the normal much tighter and
higher: max(GEV(0, 1, 0.2), N(2.5, 0.3)). Physically this is the sharper version of the same
story &mdash; ordinary snowmelt maxima cluster tightly because they are set by a fairly repeatable
snowpack, and the heavy tail is entirely rain-on-snow. Statistically it drags a global fit hard
(the MLE converges to &sigma; = 0.39 against a true 1, and is 49% low at the 1000-year level) while
leaving every level above the median GEV-governed, so a tail-weighted criterion gets a clean and
precise read: sd(&xi;&#770;) = 0.099 for the expectile version.</p>
</div>
FIG_SHARP
<div class="col">
<p>Here the composite expectile estimator wins, and the composite quantile estimator does not:</p>
</div>
<div class="tablewrap">
<table>
<caption>max(GEV(0, 1, 0.2), N(2.5, 0.3)), n = 100, 2000 replicates, against the GEV MLE.
Mean squared error is reported with its paired Monte Carlo standard error; because a handful of
extreme replicates dominate it, the median absolute error and the head-to-head win rate on the
same datasets are the more informative summaries.</caption>
<thead><tr>
<th rowspan="2">estimator</th>
<th colspan="2">MSE ratio</th><th colspan="2">median |error| ratio</th><th colspan="2">closer than the MLE</th></tr>
<tr><th>T = 200</th><th>T = 1000</th><th>T = 200</th><th>T = 1000</th><th>T = 200</th><th>T = 1000</th></tr></thead>
<tbody>
<tr><td>Composite expectile (L2), p&#8320; = 0.5</td>
  <td class="tie">0.81 <span class="se">&plusmn; 0.15</span></td><td class="tie">1.10 <span class="se">&plusmn; 0.28</span></td>
  <td class="win">0.59</td><td class="win">0.58</td><td class="win">88%</td><td class="win">84%</td></tr>
<tr><td>Composite quantile (L1), p&#8320; = 0.5</td>
  <td class="lose">1.54 <span class="se">&plusmn; 0.07</span></td><td class="lose">12.78 <span class="se">&plusmn; 0.50</span></td>
  <td>0.75</td><td class="lose">1.79</td><td>68%</td><td class="lose">38%</td></tr>
<tr><td>Composite quantile (L1), p&#8320; = 0.9</td>
  <td class="win">0.84 <span class="se">&plusmn; 0.05</span></td><td class="lose">2.82 <span class="se">&plusmn; 0.26</span></td>
  <td class="win">0.60</td><td>0.81</td><td class="win">82%</td><td>66%</td></tr>
<tr><td>GEV L-moments</td>
  <td class="win">0.72 <span class="se">&plusmn; 0.02</span></td><td class="lose">1.26 <span class="se">&plusmn; 0.11</span></td>
  <td class="win">0.70</td><td>0.76</td><td class="win">86%</td><td>77%</td></tr>
<tr class="ref"><td>GEV maximum likelihood</td><td>1.00</td><td>1.00</td><td>1.00</td><td>1.00</td><td>&mdash;</td><td>&mdash;</td></tr>
</tbody></table>
</div>
<div class="col">
<p>On typical error the composite expectile estimator is <strong>40% closer to the truth than the
MLE at every return period from 100 to 1000</strong>, and it is the closer of the two on
<strong>84% of datasets</strong> at T = 1000. The composite quantile estimator does not share
this at all &mdash; at T = 1000 it is worse than the MLE on the same criterion, winning only 38%
of the time. So this setting cleanly separates the two losses, in L2's favour, in exactly the way
the original hypothesis predicted.</p>

<p>Mean squared error is the one criterion that still calls it a tie at T = 1000
(1.10 &plusmn; 0.28). That is not a contradiction: the expectile estimator's error distribution has
a small proportion of very poor replicates, and squared error is dominated by them. Which summary
matters depends on whether an occasional large overestimate of a design flood is worse than a
systematic underestimate of every one &mdash; and in flood frequency it usually is not.</p>
</div>
</section>
"""

S7E = """
<section id="elastile">
<h2><span class="num">08</span><span>Mixing the two losses: the &alpha;-elastile</span></h2>
<div class="col">
<p>The two losses fail in opposite directions &mdash; L1 unbiased and wild, L2 steadier and
slightly biased &mdash; which invites an elastic-net-style compromise. Taking a convex combination
of the two check functions at the same level,</p>
</div>
<div class="eq">&rho;<sup>&alpha;</sup><sub>p</sub>(y, t) = (&alpha;/s) |p &minus; I(y &lt; t)| (y &minus; t)<sup>2</sup> + (1 &minus; &alpha;) (p &minus; I(y &lt; t)) (y &minus; t)</div>
<div class="col">
<p>gives a strictly consistent scoring function for a new functional &mdash; the
<strong>&alpha;-elastile</strong>. Differentiating the expected loss identifies it as the unique
root of</p>
</div>
<div class="eq">G(t) = (2&alpha;/s) [ (1 &minus; 2p) &phi;(t) + (1 &minus; p)(t &minus; m) ] + (1 &minus; &alpha;) [ F(t) &minus; p ] = 0</div>
<div class="col">
<p>The first bracket vanishes at the expectile, the second at the quantile, and both are increasing
in t &mdash; so G is strictly increasing, the root is unique, and it always lies <em>between</em>
the quantile and the expectile, which makes for a perfect bracket in the solver.
&alpha; = 0 recovers the quantile and &alpha; = 1 the expectile.</p>

<p>The scale s carries units of y and is what makes &alpha; = 1/2 mean anything: without it the L2
term is in units of y&sup2; and the L1 term in units of y. It is fixed per dataset at the ratio of
the two integrated losses evaluated at a preliminary L-moment fit, so the two contribute equally
there. Being constant during the optimisation, it does not disturb the properness of the
criterion; it only chooses the units in which &alpha; interpolates.</p>

<h3>The contamination interpolates smoothly</h3>
<p>The reason to expect the mixture to help is that the structural obstacle of section 3 arrives
gradually. The percentage by which the truth's functional exceeds its GEV component's &mdash; the
contamination no weight can remove &mdash; runs smoothly from the quantile's zero to the
expectile's:</p>
</div>
<div class="tablewrap">
<table>
<caption>Contamination surviving in the &alpha;-elastile of the contaminated truth, and the
asymptotic return-level bias each &alpha; converges to.</caption>
<thead><tr><th>&alpha;</th><th>at T = 20</th><th>at T = 50</th><th>at T = 100</th><th>at T = 1000</th><th>asymptotic bias at T = 1000</th></tr></thead>
<tbody>
<tr><td>0 &nbsp;(quantile)</td><td>0.56%</td><td>0.00%</td><td>0.00%</td><td>0.00%</td><td class="win">&minus;0.08%</td></tr>
<tr><td>0.05</td><td>1.41%</td><td>0.69%</td><td>0.60%</td><td>0.43%</td><td class="win">+0.37%</td></tr>
<tr><td>0.10</td><td>2.25%</td><td>1.33%</td><td>1.14%</td><td>0.75%</td><td class="win">+0.69%</td></tr>
<tr><td>0.35</td><td>6.21%</td><td>4.05%</td><td>3.21%</td><td>1.69%</td><td>+1.36%</td></tr>
<tr><td>0.50</td><td>8.48%</td><td>5.45%</td><td>4.20%</td><td>2.05%</td><td>+1.39%</td></tr>
<tr><td>1 &nbsp;(expectile)</td><td>15.99%</td><td>9.73%</td><td>7.03%</td><td>2.89%</td><td>+0.70%</td></tr>
</tbody></table>
</div>
<div class="col">
<p>So a small &alpha; keeps almost all of L1's asymptotic unbiasedness. The question is whether it
also picks up L2's variance reduction, and it does &mdash; monotonically. At the good weight
(p&#8320; = 0.5) and T = 1000, the standard deviation of the estimated return level falls from
13.51 at &alpha; = 0 to 6.49 at &alpha; = 1 without interruption, while the bias crosses zero
around &alpha; = 0.2. That is exactly the configuration in which an interior optimum exists.</p>
</div>
FIG_EL
<div class="tablewrap">
<table>
<caption>&alpha;-sweep at p&#8320; = 0.5, n = 100, 2000 replicates. MSE relative to the GEV MLE.
The best &alpha; is interior at every return period up to 500, and rises with T.</caption>
<thead><tr><th>&alpha;</th><th>T = 20</th><th>T = 50</th><th>T = 100</th><th>T = 200</th><th>T = 500</th><th>T = 1000</th></tr></thead>
<tbody>
<tr><td>0 &nbsp;(pure L1)</td><td>0.98</td><td>1.43</td><td>1.44</td><td>1.76</td><td>2.99</td><td>4.64</td></tr>
<tr><td>0.10</td><td>0.85</td><td>1.16</td><td>1.00</td><td>1.04</td><td>1.38</td><td>1.84</td></tr>
<tr><td>0.20</td><td class="win">0.84</td><td>1.11</td><td>0.93</td><td>0.93</td><td>1.15</td><td>1.48</td></tr>
<tr><td>0.35</td><td>0.87</td><td class="win">1.10</td><td class="win">0.90</td><td>0.88</td><td>1.02</td><td>1.26</td></tr>
<tr><td>0.50</td><td>0.96</td><td>1.12</td><td>0.90</td><td class="win">0.86</td><td>0.97</td><td>1.16</td></tr>
<tr><td>0.70</td><td>1.21</td><td>1.20</td><td>0.94</td><td>0.88</td><td class="win">0.95</td><td>1.09</td></tr>
<tr><td>1 &nbsp;(pure L2)</td><td>1.95</td><td>1.48</td><td>1.07</td><td>0.98</td><td>0.98</td><td class="win">1.05</td></tr>
<tr class="ref"><td>GEV L-moments</td><td>1.35</td><td>1.30</td><td>1.14</td><td>1.10</td><td>1.09</td><td>1.11</td></tr>
<tr class="ref"><td>gain over the better pure loss</td><td>0.85</td><td>0.77</td><td>0.84</td><td>0.89</td><td>0.97</td><td>1.00</td></tr>
</tbody></table>
</div>
<div class="col">
<p><strong>Mixing helps, by 5 to 23% over whichever pure loss is better, and the optimal &alpha;
rises with the return period</strong> &mdash; from about 0.1 at T = 10 to pure L2 at T = 1000.
That is the elastic-net pattern, and it has a clean reading here: the further out you extrapolate,
the more the variance term dominates and the more of the expectile's smoothing you want.</p>

<p>Two caveats. The benefit depends on the weight being well placed: repeating the sweep at
p&#8320; = 0.95 gives an interior optimum only at T &le; 20, and the pure losses win everywhere
beyond. And &alpha; is a tuning parameter, so a fair account of it needs the selection cost
folded in &mdash; these are oracle values.</p>
</div>
</section>
"""



SG = """
<section id="graft">
<h2><span class="num">09</span><span>Grafting the composite fit onto an empirical body</span></h2>
<div class="col">
<p>Everything above reports the composite fit as if it were the whole distribution, which is
unfair to it and unlike how it would be used. The estimator buys its tail by abandoning the body,
and the damage is severe: at p&#8320; = 0.95 the mean squared error of the 2-year return level is
<strong>eleven thousand times</strong> the MLE's for the quantile version and a thousand times for
the expectile version. That is not a subtlety, it is a ruined distribution.</p>

<p>The natural repair is to use the composite fit only as the <em>tail</em> of a smooth graft whose
body is the empirical distribution. The construction used here is the hazard-mixture smooth graft
(Coia and De Michele), taken from <code>probaverse/distplyr</code> rather than reimplemented, which
mixes body and tail on the hazard rather than the density and so has no jump or kink at the
handover and no normalising constant. Writing w for the handover weight and
<span class="mono">M = (1-w) S<sub>body</sub> + w S<sub>tail</sub></span>, the survival is
<span class="mono">S = M C</span> with a correction factor C driven by w'.</p>

<p>With a weight that is exactly 0 below x<sub>lo</sub> and exactly 1 above x<sub>hi</sub>, the
construction collapses to three regimes, which is worth spelling out because it explains the
results:</p>
<ul>
<li>below x<sub>lo</sub>: w' = 0 so C = 1, and the graft <em>is</em> the empirical distribution;</li>
<li>between: the transition, the only part needing the integral;</li>
<li>above x<sub>hi</sub>: w' = 0 again, so C freezes at a constant c and
<span class="mono">S(x) = c &middot; S<sub>tail</sub>(x)</span> &mdash; the fitted tail,
<em>rescaled by a single constant</em> that reconciles it with the body's mass at the handover.</li>
</ul>
<p>That last point matters: the graft does not merely paste the fitted tail on, it re-anchors it.
So one should expect the body to be fixed by construction, and the tail to change too.</p>
</div>
FIG_SG
<div class="col">
<p>Both expectations hold. The handover weight was run two ways &mdash; the same weight used in the
fitting criterion, as proposed, and a fixed [0.90, 0.98] window &mdash; because the two weights do
different jobs: one says which levels to fit, the other says where to stop trusting the data.</p>
</div>
<div class="tablewrap">
<table>
<caption>n = 100, 2000 replicates, MSE relative to the GEV MLE with paired Monte Carlo standard
errors. No fits failed. The last two rows are the control: grafting a <em>global</em> fit onto the
same body, to test whether any of the gain is simply the empirical body.</caption>
<thead><tr><th>fit used as the tail</th><th>T = 2</th><th>T = 10</th><th>T = 100</th><th>T = 200</th><th>T = 1000</th><th>closer than MLE at T = 1000</th></tr></thead>
<tbody>
<tr><td>Composite expectile p&#8320;=0.95, <em>ungrafted</em></td>
  <td class="lose">1168 <span class="se">&plusmn; 61</span></td><td class="lose">12.8</td><td class="lose">3.58</td><td class="lose">2.07</td><td>1.10 <span class="se">&plusmn; 0.04</span></td><td>63%</td></tr>
<tr><td>&nbsp;&nbsp;&rarr; grafted, handover 0.95&ndash;0.98</td>
  <td class="win">1.15 <span class="se">&plusmn; 0.03</span></td><td class="win">0.47</td><td>1.52</td><td>1.24</td><td class="tie">1.01 <span class="se">&plusmn; 0.02</span></td><td>59%</td></tr>
<tr><td>Composite expectile p&#8320;=0.50, <em>ungrafted</em></td>
  <td class="lose">2.98 <span class="se">&plusmn; 0.08</span></td><td>1.17</td><td>1.07</td><td>0.98</td><td>1.05 <span class="se">&plusmn; 0.05</span></td><td>78%</td></tr>
<tr><td>&nbsp;&nbsp;&rarr; grafted, handover 0.90&ndash;0.98</td>
  <td class="win">1.15 <span class="se">&plusmn; 0.03</span></td><td class="win">0.47</td><td class="win">0.92 <span class="se">&plusmn; 0.02</span></td><td class="win">0.82 <span class="se">&plusmn; 0.02</span></td><td class="win">0.84 <span class="se">&plusmn; 0.03</span></td><td class="win">80%</td></tr>
<tr><td>&nbsp;&nbsp;&rarr; grafted, handover 0.50&ndash;0.53 (same as fitting weight)</td>
  <td class="win">1.15</td><td class="lose">1.89</td><td>1.32</td><td>1.20</td><td>1.34 <span class="se">&plusmn; 0.08</span></td><td>76%</td></tr>
<tr><td>Composite quantile p&#8320;=0.50 &rarr; grafted, 0.90&ndash;0.98</td>
  <td class="win">1.15</td><td class="win">0.47</td><td class="lose">2.01</td><td class="lose">2.46</td><td class="lose">6.38 <span class="se">&plusmn; 0.58</span></td><td>58%</td></tr>
<tr class="ref"><td>Empirical distribution alone</td>
  <td>1.15</td><td>0.47</td><td>5.37</td><td>2.03</td><td>1.13</td><td>43%</td></tr>
<tr class="ref"><td>GEV MLE &rarr; grafted, 0.90&ndash;0.98 <em>(control)</em></td>
  <td>1.15</td><td>0.47</td><td>1.13 <span class="se">&plusmn; 0.01</span></td><td>1.07 <span class="se">&plusmn; 0.01</span></td><td>1.02 <span class="se">&plusmn; 0.01</span></td><td>50%</td></tr>
<tr class="ref"><td>GEV L-moments &rarr; grafted, 0.90&ndash;0.98 <em>(control)</em></td>
  <td>1.15</td><td>0.47</td><td>1.24</td><td>1.15</td><td>1.14</td><td>62%</td></tr>
</tbody></table>
</div>
<div class="col">
<p><strong>The body is completely repaired.</strong> Eleven thousand becomes 1.15, which is the
empirical distribution's own figure &mdash; as it must be, since below the handover the graft is
exactly the empirical distribution. The composite estimator's worst liability simply disappears.</p>

<p><strong>The tail improves as well</strong>, which was not guaranteed. The best combination,
composite expectile fitted at p&#8320; = 0.5 and handed over on [0.90, 0.98], beats the GEV MLE at
every return period from 50 to 1000 &mdash; 0.92, 0.82, 0.80, 0.84 &mdash; each five to eight
standard errors below one, and it is closer to the truth than the MLE on <strong>80% of
datasets</strong>. Ungrafted, the same fit was 1.07, 0.98, 0.98, 1.05: a tie. The re-anchoring
constant c is doing real work.</p>

<div class="note"><span class="lab">The control settles it</span>
<p>Grafting the <em>MLE</em> onto the same empirical body gives 1.13, 1.07, 1.03, 1.02 and a 50%
win rate &mdash; that is, nothing. Grafting L-moments is likewise neutral. So the gain is not the
empirical body, and it is not the graft; it is the composite expectile tail, which only becomes
usable once the graft takes the body off its hands. Neither ingredient works alone.</p></div>

<p>Two practical points. The handover weight should <em>not</em> be tied to the fitting weight:
handing over at [0.50, 0.53] because the criterion was fitted from p&#8320; = 0.5 means trusting
the parametric model from the 53rd percentile up, and costs a factor of 1.6 at T = 1000 against
the decoupled [0.90, 0.98] version. And the graft does not rescue the L1 fit &mdash; grafted, the
composite quantile estimator is still 6.4&times; the MLE at T = 1000. The tail model has to be
worth grafting.</p>
</div>
</section>
"""


if __name__ == "__main__":
    html = build()
    (ROOT / "report" / "report.html").write_text(html)
    print("wrote report/report.html  (%.1f MB)" % (len(html) / 1e6))

S10 = """
<section id="anchoring-fix">
<h2><span class="num">10</span><span>Fixing the mean anchoring: it works, and it does not pay</span></h2>
<div class="col">
<p>Section 3 identified the one structural disadvantage L2 carries against L1: the expectile is
anchored to the distribution's mean, so body contamination never fully leaves it. The
identification equation splits cleanly, though &mdash;</p>
</div>
<div class="eq">k &phi;(x) + m &minus; x = 0</div>
<div class="col">
<p>&mdash; into a tail-local part (the partial moment &phi;, where a GEV fitted to a contaminated
sample is right) and a global part (the mean, where it is not). So take the second from the data:
replace the model's mean by the <em>sample</em> mean and leave &phi; as the model's. The mean
becomes a nuisance parameter profiled out nonparametrically.</p>

<p>As a piece of statistics it works exactly as designed. The contamination surviving in the
expectile at p = 0.98 falls from 9.73% to 0.0001%, and the estimator's asymptotic target moves onto
the true tail parameters:</p>
</div>
<div class="tablewrap">
<table>
<caption>Asymptotic targets, fitted to a sample of 2&times;10<sup>6</sup>. The true tail is
(0, 1, 0.20).</caption>
<thead><tr><th>estimator</th><th>&mu;</th><th>&sigma;</th><th>&xi;</th><th>bias at T = 100</th><th>at T = 1000</th></tr></thead>
<tbody>
<tr><td>Composite expectile, p&#8320; = 0.95</td><td>0.836</td><td>0.872</td><td>0.220</td><td>+2.86%</td><td>+0.51%</td></tr>
<tr><td>&nbsp;&nbsp;&rarr; mean-anchored</td><td class="win">0.026</td><td class="win">0.989</td><td class="win">0.202</td><td class="win">&minus;0.32%</td><td class="win">&minus;0.23%</td></tr>
<tr><td>Composite expectile, p&#8320; = 0.90</td><td>1.006</td><td>0.797</td><td>0.237</td><td>+1.22%</td><td>+0.25%</td></tr>
<tr><td>&nbsp;&nbsp;&rarr; mean-anchored</td><td class="win">0.161</td><td class="win">0.953</td><td class="win">0.208</td><td class="win">&minus;0.49%</td><td class="win">&minus;0.26%</td></tr>
</tbody></table>
</div>
<div class="col">
<p>The structural obstacle is gone: anchored, the L2 estimator is as asymptotically unbiased as L1.</p>

<div class="note"><span class="lab">And yet it is worse</span>
<p>At n = 100 the correction loses on every criterion at every return period. The median fitted
shape goes from 0.209 (ordinary, p&#8320; = 0.5) to 0.079, and at p&#8320; = 0.95 it collapses to
&minus;0.45, the shape bound. Mean squared error at T = 1000 goes from 1.05 to 1.09 relative to the
MLE; median absolute error from 4.68 to 4.81; the head-to-head win rate from 78% to 76%. Grafted,
the ordinary version still wins (0.84&times; the MLE against 0.92&times;).</p></div>

<p>The reason is visible in the equation. An error &delta; in the anchor moves the solution by
about &delta;/(1 + k S(x)), which for large p is essentially &delta; itself &mdash; the whole
expectile curve translates by the sample mean's estimation error. Under a contaminated,
heavy-tailed truth at n = 100 that error is not small, and it is injected at every level. Worse,
the criterion must now reconcile a data-driven anchor with a model-driven &phi;, and when the
weight is confined to the tail the model has little left to adjust except the shape &mdash; which
is why &xi;&#770; is driven to its bound.</p>

<p>So this is a clean negative result, and it closes off the direction rather than opening it: the
mean anchoring is real, the fix removes it exactly, and the fix is not worth its variance. If the
idea is to be revived it needs a <em>low-variance</em> estimate of the anchor &mdash; a shrunk or
model-averaged mean, or the mean of a trimmed body &mdash; not the raw sample mean.</p>
</div>
</section>
"""

S11 = """
<section id="handoff">
<h2><span class="num">11</span><span>Where the handover should sit</span></h2>
<div class="col">
<p>Section 9 used a fixed [0.90, 0.98] handover on the grounds that the fitting weight and the
handover weight do different jobs. Sweeping it further out tests that: if the empirical body is
good, handing over later should be better.</p>

<p>It is not. Later is monotonically worse:</p>
</div>
<div class="tablewrap">
<table>
<caption>Composite expectile fitted at p&#8320; = 0.5, grafted onto an empirical body with the
handover swept. MSE relative to the GEV MLE, n = 100, 2000 replicates, paired standard errors.</caption>
<thead><tr><th>handover</th><th>T = 20</th><th>T = 50</th><th>T = 100</th><th>T = 200</th><th>T = 500</th><th>T = 1000</th></tr></thead>
<tbody>
<tr><td>0.90 &ndash; 0.98</td><td class="win">1.47</td><td class="win">1.28</td><td class="win">0.92 <span class="se">&plusmn;0.02</span></td><td class="win">0.82 <span class="se">&plusmn;0.02</span></td><td class="win">0.80 <span class="se">&plusmn;0.03</span></td><td class="win">0.84 <span class="se">&plusmn;0.03</span></td></tr>
<tr><td>0.95 &ndash; 0.99</td><td>1.79</td><td>1.71</td><td>1.16</td><td>1.00</td><td>0.95</td><td>1.00 <span class="se">&plusmn;0.04</span></td></tr>
<tr><td>0.97 &ndash; 0.995</td><td>1.79</td><td>1.94</td><td class="lose">1.43</td><td class="lose">1.81</td><td class="lose">1.83</td><td class="lose">2.00 <span class="se">&plusmn;0.13</span></td></tr>
<tr><td>0.99 &ndash; 0.999</td><td>1.79</td><td>1.93</td><td class="lose">1.81</td><td class="lose">2.33</td><td class="lose">2.32</td><td class="lose">2.51 <span class="se">&plusmn;0.18</span></td></tr>
</tbody></table>
</div>
<div class="col">
<p>The mechanism is the empirical distribution's own limit. Its quantile function saturates at the
sample maximum, so it carries no information above about p = 1 &minus; 1/n &mdash; 0.99 at n = 100.
A handover that only <em>starts</em> at 0.99 leaves the graft relying on the empirical body across
levels the empirical body cannot see, and the return-level curve flattens onto order statistics
near the maximum long before the parametric tail is allowed to take over.</p>

<p>That gives a concrete design rule, and it is the mirror image of the one in section 7. The
weight must not extend past where the data can speak; and the <em>handover</em> must complete
before the empirical body runs out. With n observations, both point to the same place: finish the
handover by roughly 1 &minus; 2/n. At n = 100 that is 0.98, which is exactly where the best
setting lands.</p>

<p>So the two weights should be decoupled, but not in the direction of pushing the handover later.
The fitting weight can start anywhere the model is trustworthy; the handover is pinned near the
top of the data by the sample size, not by the fit.</p>
</div>
</section>
"""

S12 = """
<section id="alphap">
<h2><span class="num">12</span><span>Letting &alpha; depend on the level</span></h2>
<div class="col">
<p>Since the best &alpha; rises with the return period, the obvious next move is to let &alpha;
depend on p. Does that still elicit anything?</p>

<p><strong>Yes, and the argument is short.</strong> Properness here is a pointwise-in-p property.
For each fixed level p, the mixed check function &rho;<sup>&alpha;</sup><sub>p</sub> is strictly
consistent for the &alpha;-elastile at that level, because the expected loss is strictly convex in
the prediction with a unique minimiser. Nothing in that argument cares whether the &alpha; used at
level p is the same as the one used at level p&prime;. So for any measurable schedule
&alpha;(&middot;), the criterion</p>
</div>
<div class="eq">R(&theta;) = &#8747; w(p) E[ &rho;<sup>&alpha;(p)</sup><sub>p</sub>(Y, T(p|&theta;)) ] dp</div>
<div class="col">
<p>is bounded below by the sum of the pointwise minima, with equality exactly when
T(p|&theta;) = T<sub>&alpha;(p)</sub>(p; F) for w-almost every p. It is a strictly consistent
scoring rule for the <em>curve</em> p &#8614; T<sub>&alpha;(p)</sub>(p; F). Fisher consistency
survives too: under correct specification over the weighted region, the truth's curve <em>is</em>
the model's curve, for any schedule.</p>

<p>Three things are worth saying about what you get.</p>

<p><strong>The target is a hybrid.</strong> Different points on the elicited curve are different
functionals of F &mdash; the 0.2-elastile at p = 0.9, the 0.8-elastile at p = 0.999. That is
perfectly well defined but it no longer has a name like "the expectile function", and any theory
you would want to import (asymptotic variance for expectile regression, say) does not transfer
without redoing it.</p>

<p><strong>Keep &alpha;(p) continuous.</strong> This is the one real constraint, and it is not
about properness &mdash; it is about attainability. A model's quantile, expectile and elastile
curves are all increasing in p. The target curve need not be. In the tail the expectile sits
<em>below</em> the quantile, so raising &alpha; pulls the target down while raising p pushes it up.
With any continuous schedule tested &mdash; &alpha; = p, a smoothstep ramp over [0.9, 0.999], even
&alpha; = 1 &minus; (1&minus;p)<sup>0.05</sup>, which is above 0.7 by p = 0.5 &mdash; the target
stayed strictly increasing, at &xi; = 0.2 and at &xi; = 0.45. A schedule that <em>jumps</em>
(&alpha; = 0 below p = 0.99, 1 above) makes the target fall by 1.27 at the jump. The criterion is
still proper there; it is simply aiming at a curve no model can match, which is misspecification
you have manufactured for yourself.</p>

<p><strong>It generalises further than &alpha;.</strong> The same pointwise argument licenses
letting <em>any</em> parameter of the level-p loss vary with p &mdash; the scale s included &mdash;
and it is also why the composite construction tolerates a different loss family at each level. The
cost is always the same: more tuning surface, and oracle-selected schedules that need their
selection cost charged before the reported gains mean anything.</p>
</div>
</section>
"""
