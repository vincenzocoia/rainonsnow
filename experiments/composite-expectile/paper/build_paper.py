"""Builds the manuscript as a single self-contained HTML file."""
import base64
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "out"


def fig(name, num, caption):
    data = base64.b64encode((OUT / name).read_bytes()).decode()
    return (f'<figure class="wide">\n'
            f'  <img src="data:image/png;base64,{data}" alt="{caption[:110]}">\n'
            f'  <figcaption><span class="fignum">Figure {num}.</span> {caption}</figcaption>\n'
            f'</figure>\n')


HEAD = """<title>Tail-Weighted Composite Losses</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=Literata:opsz,wght@7..72,400;7..72,600&display=swap">
<style>
:root{
  --ground:#fafbfc; --surface:#ffffff; --surface-2:#eef2f5;
  --ink:#141d24; --ink-2:#4b5f6d; --ink-3:#7d909c;
  --rule:#d5dee4; --rule-strong:#a9bac5;
  --accent:#1c5d7d; --accent-2:#9c4a2f; --accent-3:#2d6a4a;
  --sans:"IBM Plex Sans",ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;
  --serif:"Literata",Georgia,"Times New Roman",serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --measure:68ch; --wide:1080px;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --ground:#0e151a; --surface:#151f26; --surface-2:#1b272f;
  --ink:#e2eaf0; --ink-2:#9aaebb; --ink-3:#71858f;
  --rule:#26343d; --rule-strong:#3b4d59;
  --accent:#6bb2d6; --accent-2:#dd937f; --accent-3:#6cc79a;
}}
:root[data-theme="dark"]{
  --ground:#0e151a; --surface:#151f26; --surface-2:#1b272f;
  --ink:#e2eaf0; --ink-2:#9aaebb; --ink-3:#71858f;
  --rule:#26343d; --rule-strong:#3b4d59;
  --accent:#6bb2d6; --accent-2:#dd937f; --accent-3:#6cc79a;
}
*{box-sizing:border-box}
body{background:var(--ground);color:var(--ink);font-family:var(--serif);
  font-size:17px;line-height:1.66;-webkit-font-smoothing:antialiased}
.page{max-width:var(--wide);margin:0 auto;padding:0 28px 110px}
.col{max-width:var(--measure)}

/* --- title block ------------------------------------------------------- */
header.title{padding:64px 0 26px;border-bottom:1px solid var(--rule-strong);margin-bottom:34px}
.kicker{font-family:var(--sans);font-size:11px;font-weight:600;letter-spacing:.15em;
  text-transform:uppercase;color:var(--accent);margin-bottom:16px}
h1{font-family:var(--sans);font-weight:600;font-size:clamp(27px,3.9vw,40px);
  line-height:1.16;letter-spacing:-.02em;margin:0 0 18px;max-width:30ch;text-wrap:balance}
.sub{font-size:19px;line-height:1.55;color:var(--ink-2);max-width:var(--measure);margin:0}
.abstract{background:var(--surface);border:1px solid var(--rule);border-left:3px solid var(--accent);
  padding:22px 26px;margin:30px 0 18px;max-width:var(--measure);font-size:16px}
.abstract h2{font-family:var(--sans);font-size:11px;font-weight:600;letter-spacing:.14em;
  text-transform:uppercase;color:var(--ink-3);margin:0 0 10px;display:block}
.abstract p{margin:0 0 12px}.abstract p:last-child{margin:0}
.keywords{font-family:var(--sans);font-size:13px;color:var(--ink-2);max-width:var(--measure);
  margin:0 0 8px}
.keywords b{color:var(--ink);font-weight:600}

/* --- structure --------------------------------------------------------- */
h2{font-family:var(--sans);font-weight:600;font-size:21px;letter-spacing:-.01em;
  line-height:1.25;margin:52px 0 4px;text-wrap:balance;display:flex;gap:14px;align-items:baseline}
h2 .n{font-family:var(--mono);font-size:14px;font-weight:500;color:var(--accent);flex:none}
h3{font-family:var(--sans);font-weight:600;font-size:16px;margin:30px 0 2px;color:var(--ink)}
h3 .n{font-family:var(--mono);font-size:13px;font-weight:500;color:var(--ink-3);margin-right:9px}
p{margin:0 0 14px}h2+p,h3+p{margin-top:10px}
a{color:var(--accent);text-decoration-thickness:1px;text-underline-offset:2px}
strong{font-weight:600}
em.term{font-style:italic;color:var(--ink)}
code{font-family:var(--mono);font-size:.86em;background:var(--surface-2);padding:1px 4px;border-radius:2px}
ul,ol{max-width:var(--measure);padding-left:22px;margin:0 0 14px}li{margin-bottom:7px}

/* --- equations --------------------------------------------------------- */
.eq{display:flex;align-items:center;gap:18px;margin:20px 0;max-width:var(--measure)}
.eq .body{flex:1;font-family:var(--mono);font-size:14px;line-height:1.8;
  background:var(--surface);border-left:3px solid var(--rule-strong);
  padding:14px 18px;overflow-x:auto;white-space:pre;color:var(--ink)}
.eq .tag{font-family:var(--mono);font-size:13px;color:var(--ink-3);flex:none}
.eq .hl{color:var(--accent-2);font-weight:500}

/* --- propositions ------------------------------------------------------ */
.prop{max-width:var(--measure);background:var(--surface);border:1px solid var(--rule);
  border-top:3px solid var(--accent);padding:18px 22px;margin:24px 0}
.prop .lab{font-family:var(--sans);font-size:11px;font-weight:600;letter-spacing:.11em;
  text-transform:uppercase;color:var(--accent);display:block;margin-bottom:7px}
.prop p:last-child{margin-bottom:0}
.note{max-width:var(--measure);border-left:3px solid var(--accent-2);
  padding:4px 0 4px 18px;margin:22px 0;color:var(--ink-2);font-size:16px}
.note p:last-child{margin-bottom:0}

/* --- tables ------------------------------------------------------------ */
.tablewrap{overflow-x:auto;margin:24px 0;max-width:var(--wide)}
table{border-collapse:collapse;font-family:var(--sans);font-size:13.5px;
  font-variant-numeric:tabular-nums;min-width:100%}
caption{text-align:left;font-family:var(--sans);font-size:13px;color:var(--ink-2);
  padding-bottom:10px;max-width:var(--measure);line-height:1.5}
caption b{color:var(--ink);font-weight:600}
th,td{padding:7px 13px;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left;padding-left:0}
thead th{font-weight:600;font-size:12px;letter-spacing:.03em;color:var(--ink-2);
  border-bottom:1.5px solid var(--rule-strong)}
tbody td{border-bottom:1px solid var(--rule)}
tbody tr:last-child td{border-bottom:1.5px solid var(--rule-strong)}
tr.ref td{color:var(--ink-2);font-style:italic}
tr.sub td{color:var(--ink-3);font-size:.92em}
td.win{color:var(--accent-3);font-weight:600}
td.lose{color:var(--accent-2)}

/* --- figures ----------------------------------------------------------- */
figure{margin:30px 0;max-width:var(--wide)}
figure img{width:100%;display:block;background:#fff;border:1px solid var(--rule);border-radius:2px}
figcaption{font-family:var(--sans);font-size:13px;line-height:1.55;color:var(--ink-2);
  margin-top:10px;max-width:var(--measure)}
.fignum{color:var(--ink);font-weight:600}

/* --- references -------------------------------------------------------- */
.refs{max-width:var(--measure);font-size:15px;line-height:1.55}
.refs li{margin-bottom:11px;color:var(--ink-2)}
.refs li b{color:var(--ink);font-weight:600}
footer{margin-top:70px;padding-top:20px;border-top:1px solid var(--rule);
  font-family:var(--sans);font-size:13px;color:var(--ink-3);max-width:var(--measure)}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
@media (max-width:640px){.eq{flex-direction:column;align-items:stretch;gap:6px}
  .eq .tag{text-align:right}}
</style>

<div class="page">
<header class="title">
  <div class="kicker">Statistics of extremes &middot; manuscript draft</div>
  <h1>Composite M-quantile estimation for extreme-value models under body misspecification</h1>
  <p class="sub">A tail-weighted loss buys a real reduction in far-tail error when the body of the
  distribution is wrong, and costs 12&ndash;35% when it is not. Both numbers matter, and the
  standard comparison reports only the first.</p>
</header>
<main>

<div class="abstract">
<h2>Abstract</h2>
<p>Extreme-value models are justified asymptotically but fitted at finite block sizes and
thresholds, where the family is typically adequate in the tail and wrong in the body. Global
criteria &mdash; likelihood, L-moments &mdash; distribute that misspecification across the whole
distribution and bias the extrapolation of interest. We study estimators that integrate an
asymmetric loss over probability levels against a measure concentrated in the tail, and place them
in the M-quantile family. Replacing the pinball loss of the composite quantile estimator by its
asymmetric-least-squares analogue reduces mean squared error at the 1000-year return level by a
factor of two to six, at the cost of a structural obstacle: the expectile identification equation
contains the distribution's mean, so no weight makes the estimator asymptotically unbiased under
body contamination. A one-sided inverted-Huber influence function removes the mean exactly, giving
a functional that depends on the distribution only above a knot, and it is the most accurate
far-tail estimator we examine. Against that, we report what the approach costs when the model is
correct &mdash; 12&ndash;32% in the far tail of a correctly specified GEV, more at short return
periods, with bias in the unsafe direction &mdash; and show that 68% of the apparent gain over
peaks-over-threshold in a misspecified study is present with no misspecification at all, because
the conventional reference discards nine tenths of the sample.</p>
</div>
<p class="keywords"><b>Keywords:</b> tail robustness &middot; expectiles &middot; M-quantiles
&middot; generalized extreme value &middot; generalized Pareto &middot; return level &middot;
model misspecification &middot; flood frequency</p>
"""

S1 = """
<h2><span class="n">1</span><span>Introduction</span></h2>
<div class="col">
<p>Extreme-value practice rests on two limit theorems. Block maxima converge to the generalized
extreme value (GEV) family; exceedances over a high threshold converge to the generalized Pareto
(GPD). Both are statements about a limit, and both are applied at finite block sizes and finite
thresholds, where convergence is incomplete. The resulting misspecification is not uniform across
the distribution. In a great many applications the fitted family is a reasonable description of the
upper tail and a poor description of the body, because the body is generated by a different
physical mechanism.</p>

<p>Rain-on-snow flooding is a clean instance. Annual maxima at a snow-influenced catchment are
drawn from two populations: ordinary snowmelt freshets, which are moderate and tightly clustered,
and rain-on-snow events, which are rare and heavy. The maximum of the two has a tail governed by
the rain-driven mechanism and a body governed by the snowmelt one. A single GEV cannot describe
both, and a criterion that weights all parts of the distribution equally &mdash; maximum
likelihood, or L-moments &mdash; spends its fit reconciling them. The cost falls on the
extrapolation that the analysis exists to produce.</p>

<p>The literature offers three responses. <em class="term">Threshold selection</em> discards the
body entirely, trading bias for the variance of a smaller sample. <em class="term">Robust
estimation</em> bounds the influence of individual observations: Dupuis (1998) and Ju&aacute;rez
and Schucany (2004) for the GPD, and the minimum density power divergence estimator of Basu et al.
(1998), which is indexed by a constant running from maximum likelihood at zero to bounded influence
above it. That family shares the elastile's structure &mdash; a single dial from efficient to
robust &mdash; but runs in the opposite direction, downweighting large observations as
contamination where the tail problem treats them as signal. The third response, and the one this
paper belongs to, is <em class="term">tail-robustness</em>: weight the fitting criterion towards
the tail so that contamination in the body does not propagate into the extrapolation. Fung (2022)
proves that a weighted likelihood of this kind yields a consistent tail index under model
misspecification, motivated by insurance claim severities.</p>

<p>We study a different construction in the same spirit. Rather than weighting observations, we
integrate an asymmetric loss over <em>probability levels</em>, against a measure concentrated near
one. With the pinball loss this gives the composite quantile estimator, which is asymptotically
unbiased under body contamination but suffers a severe variance penalty in the far tail. The
asymmetric-least-squares analogue &mdash; the composite <em>expectile</em> estimator &mdash; has
not been examined, and it is the natural candidate: expectiles use the magnitudes of observations
rather than their ranks, which is where the information about the shape parameter lives.</p>

<p>The paper makes four contributions. Section 2 defines the composite M-quantile family and shows
that its members, including a convex combination of the two losses, are M-quantiles in the sense of
Breckling and Chambers (1988), with closed-form identification equations for the GEV and GPD.
Section 3 characterises when the criterion is well posed: which weights are admissible, and which
moment conditions each loss requires. Section 4 establishes the central obstacle &mdash; the
expectile is anchored to the distribution's mean, so no weight removes body contamination
asymptotically &mdash; and gives an influence function that removes it structurally. Section 5
reports a simulation study designed to separate what the criterion contributes from what the
comparison contributes, and finds that most of the apparent advantage over conventional practice is
the latter.</p>
</div>
"""

S2 = """
<h2><span class="n">2</span><span>The composite M-quantile family</span></h2>
<h3><span class="n">2.1</span>The criterion</h3>
<div class="col">
<p>Let <code>T(p|&theta;)</code> denote a level-p functional of a parametric family indexed by
&theta; &mdash; its quantile function, expectile function, or another member of the family below.
Let w be a weight on (0, 1), non-negative and integrable. Given observations
y&#8321;, &hellip;, y&#8345; we minimise</p>
</div>
<div class="eq"><div class="body">S(&theta;) = &Sigma;&#8202;<sub>i</sub> &#8747;&#8320;&sup1;  w(p) <span class="hl">|</span>p &minus; I(y&#8202;<sub>i</sub> &lt; T(p|&theta;))<span class="hl">|</span> &rho;(y&#8202;<sub>i</sub> &minus; T(p|&theta;)) dp</div><div class="tag">(1)</div></div>
<div class="col">
<p>with &rho;(u) = |u| recovering the composite quantile estimator and &rho;(u) = u&sup2; giving the
composite expectile estimator. The absolute value on the asymmetry factor is essential and is
easily dropped: without it the factor (p &minus; I) is negative for every observation below
T(p|&theta;), so the criterion is unbounded below and a numerical optimiser runs to whichever
boundary it is given. The pinball loss (p &minus; I)(y &minus; q) carries no absolute value because
it is already non-negative; the asymmetric square does not share that property.</p>

<p>With the absolute value in place, each integrand is the canonical strictly consistent scoring
function for its functional, so the population criterion is minimised by a &theta; whose level-p
functional matches the truth's wherever w &gt; 0. The construction therefore inherits elicitability
pointwise in p, and any weight yields a proper criterion.</p>
</div>

<h3><span class="n">2.2</span>Members</h3>
<div class="col">
<p>All of these estimators are M-estimators, with
&rho;&#771;(y;&theta;) = &#8747; w(p) &rho;<sub>p</sub>(y, T(p|&theta;)) dp, so standard sandwich
asymptotics apply without modification. More usefully, their <em>functionals</em> sit in one family.
The level-p M-quantile of Breckling and Chambers (1988) for an influence function &psi; is the
root of</p>
</div>
<div class="eq"><div class="body">E[ |p &minus; I(Y &lt; t)| &psi;(Y &minus; t) ] = 0</div><div class="tag">(2)</div></div>
<div class="col">
<p>with &psi;(u) = sign(u) giving the quantile and &psi;(u) = u the expectile. A convex combination
of the two check functions &mdash; which we call the <em class="term">&alpha;-elastile</em>, with
scale s chosen to put the two losses on comparable footing &mdash; is the M-quantile with</p>
</div>
<div class="eq"><div class="body">&psi;(u) = (2&alpha;/s)&#8202;u + (1 &minus; &alpha;)&#8202;sign(u)</div><div class="tag">(3)</div></div>
<div class="col">
<p>This is not Huber's influence function, which switches between the two at a knot; here they are
added everywhere. Figure 1 shows the family. It also shows why the classical literature never
pursued the additive mixture: &psi; in (3) is unbounded for every &alpha; &gt; 0, so the estimator
has infinite gross-error sensitivity and zero breakdown point. Judged by the criterion that
organises robust statistics it fails at the first screen. That criterion is the wrong one here. In
tail estimation the large observations are not contamination; they carry nearly all the information
about the shape parameter, and bounding their influence discards it deliberately.</p>
</div>
FIG_PSI

<h3><span class="n">2.3</span>Closed forms</h3>
<div class="col">
<p>Equation (2) is tractable for both extreme-value families because each identification equation
reduces to the first partial moment</p>
</div>
<div class="eq"><div class="body">&phi;(x) = E[(Y &minus; x)<sup>+</sup>]</div><div class="tag">(4)</div></div>
<div class="col">
<p>which is available in closed form: for the GPD,
&phi;(x) = (&sigma; + &xi;(x &minus; &mu;))&#8202;S(x)/(1 &minus; &xi;); for the GEV, in terms of a
lower incomplete gamma function. The expectile is then the unique root of</p>
</div>
<div class="eq"><div class="body">k&#8202;&phi;(t) + m &minus; t = 0,    k = (2p &minus; 1)/(1 &minus; p),   m = E[Y]</div><div class="tag">(5)</div></div>
<div class="col">
<p>a strictly decreasing function of t admitting a safeguarded Newton solution in a handful of
iterations. The elastile requires only the same ingredients. In consequence the composite criterion
can be evaluated on a Gauss&ndash;Legendre grid over levels without any nested numerical
integration, and the sums over observations reduce to partial sums of y and y&sup2;, so the whole
grid costs one pass over two cumulative sums.</p>
</div>
"""

S3 = """
<h2><span class="n">3</span><span>When the criterion is well posed</span></h2>
<h3><span class="n">3.1</span>Admissible weights</h3>
<div class="col">
<p>The weight cannot be arbitrary. For a weight diverging polynomially at the upper end,
w(p) &prop; (1 &minus; p)<sup>&minus;a</sup>, the population criterion is finite if and only if
a &lt; 2 &minus; &xi; in the L1 case and a &lt; 2 &minus; 2&xi; in the L2 case. The tail of the
integrand is governed by the growth of the functional, which is (1 &minus; p)<sup>&minus;&xi;</sup>
for the quantile and the same order for the expectile; squaring the residual doubles the exponent.
A weight decaying too slowly makes the criterion infinite at every &theta;, and the estimator is
then not merely inefficient but undefined.</p>
</div>

<h3><span class="n">3.2</span>Moment conditions</h3>
<div class="col">
<p>A second requirement is on the truth rather than the weight, and it separates into three
questions that are easily conflated. Writing the loss as
&rho;<sub>p,a</sub>(y, t) = |p &minus; I(y &lt; t)|&#8202;|y &minus; t|<sup>a</sup>: the raw
criterion needs E|Y|<sup>a</sup>; the loss <em>difference</em> against a fixed reference function,
which has the same minimiser, is bounded by C(1 + |Y|<sup>a&minus;1</sup>), as is the estimating
equation it differentiates to, so the functional is defined under E|Y|<sup>a&minus;1</sup>; and the
sandwich variance needs the square of that score, E|Y|<sup>2(a&minus;1)</sup>. Reading these off a
GPD tail with index 1/&xi; gives Table 1.</p>
</div>
<div class="tablewrap">
<table>
<caption><b>Table 1.</b> What each loss requires of the tail, for a GPD with shape &xi;. Only the
middle column decides whether the target functional exists; the first is repaired for free by the
loss-difference form.</caption>
<thead><tr><th>loss</th><th>raw criterion finite</th><th>functional defined</th><th>root-n asymptotics</th></tr></thead>
<tbody>
<tr class="ref"><td>maximum likelihood</td><td>&mdash;</td><td>any &xi;</td><td>&xi; &gt; &minus;1/2</td></tr>
<tr class="ref"><td>L-moments</td><td>&xi; &lt; 1</td><td>&xi; &lt; 1</td><td class="lose">&xi; &lt; 1/2</td></tr>
<tr><td>pinball (a = 1)</td><td>&xi; &lt; 1</td><td class="win">any &xi;</td><td class="win">any &xi;</td></tr>
<tr><td>asymmetric square (a = 2)</td><td class="lose">&xi; &lt; 1/2</td><td>&xi; &lt; 1</td><td class="lose">&xi; &lt; 1/2</td></tr>
<tr><td>&alpha;-elastile, any &alpha; &gt; 0</td><td class="lose">&xi; &lt; 1/2</td><td>&xi; &lt; 1</td><td class="lose">&xi; &lt; 1/2</td></tr>
</tbody></table>
</div>
<div class="col">
<p>Two consequences. The quantile requires nothing: the &xi; &lt; 1 in its first column is an
artefact of not recentring. And the elastile buys no additional domain &mdash; a sum is finite only
if both terms are, so every &alpha; above zero inherits the expectile's thresholds. Mixing trades
bias against variance <em>within</em> a domain and cannot widen it. Note also that L-moments and
the asymmetric square have identical requirements, so the practitioner comparison in Section 5 is
like-for-like on tail heaviness.</p>
</div>

<h3><span class="n">3.3</span>Weighting beyond the data</h3>
<div class="col">
<p>A third and purely finite-sample trap: a weight left at one up to p = 1 places a substantial
share of its mass at levels above the largest observation. There the empirical functional has
saturated at the sample maximum, so the criterion can only pull the fitted tail down. With a weight
identically one above the 95th percentile, 45% of it sits above the sample maximum at n = 50 and
27% at n = 100, and the median fitted shape collapses to &minus;0.82 in the former case against a
true 0.2. The fitted shape
tracks that share almost exactly, and both converge to the truth as n grows, so this is a
finite-sample defect of the weight rather than of the estimator. It nonetheless dictates practice:
the weight must be chosen relative to the sample size, not in the abstract.</p>
</div>
"""

S4 = """
<h2><span class="n">4</span><span>Mean anchoring, and a loss that removes it</span></h2>
<div class="col">
<p>The expectile's advantage over the quantile &mdash; that it uses magnitudes &mdash; is also the
source of a structural problem, visible in (5).</p>
</div>
<div class="prop"><span class="lab">Proposition 1 &mdash; mean anchoring</span>
<p>The level-p expectile solves p&#8202;E[(Y &minus; t)<sup>+</sup>] =
(1 &minus; p)&#8202;E[(t &minus; Y)<sup>+</sup>], and the lower partial moment satisfies
E[(t &minus; Y)<sup>+</sup>] = t &minus; E[Y] + E[(Y &minus; t)<sup>+</sup>]. The functional
therefore depends on the distribution through E[Y] as well as through its tail. Let F&#8320; be a
model and F a contaminated version agreeing with F&#8320; above some x&#8320; but differing below.
Then the expectiles of F and F&#8320; differ at <em>every</em> level p &lt; 1, however far into the
tail, because E[Y] differs. No choice of weight w in (1) removes this: the criterion is minimised
by matching a functional that is itself contaminated.</p></div>
<div class="col">
<p>The quantile has no such term &mdash; it depends on F only through F(t) &mdash; which is why the
composite quantile estimator is asymptotically unbiased under body contamination and the composite
expectile estimator is not. The effect is not negligible. For the contaminated GEV of Section 5,
the truth's expectile at p = 0.98 sits 9.73% above its own GEV component's, and at p = 0.99 still
7.03%, while the corresponding quantile figures are zero to four decimals.</p>

<p>One obvious repair is to substitute the sample mean for the model's in (5). This works
asymptotically &mdash; it reduces the surviving contamination at p = 0.98 from 9.73% to 0.0001% and
puts the estimator's asymptotic target on the true parameters &mdash; and fails in practice, because
at n = 100 the sample mean of a heavy-tailed sample is too noisy: the corrected estimator loses on
every finite-sample criterion we measured. A structural repair is available instead.</p>
</div>
<div class="prop"><span class="lab">Proposition 2 &mdash; a functional of the upper tail alone</span>
<p>Take the one-sided influence function &psi;<sub>c</sub>(u) = max(u, &minus;c), linear above the
fitted level and capped at &minus;c below it. Substituting into (2) and using the identity above,
the two mean terms cancel and the identification equation becomes</p>
<div class="eq" style="margin:14px 0"><div class="body">p&#8202;&phi;(t) = (1 &minus; p)&#8202;[ c + &phi;(t) &minus; &phi;(t &minus; c) ]</div><div class="tag">(6)</div></div>
<p>which involves F only through &phi; at t and at t &minus; c. Since
&phi;(x) = &#8747;<sub>x</sub><sup>&infin;</sup> S(u) du, the level-p functional
t<sub>p</sub>(F) depends on F only through its restriction to [t<sub>p</sub> &minus; c, &infin;). It
is therefore invariant to any rearrangement of probability mass strictly below
t<sub>p</sub> &minus; c that preserves the total mass there.</p></div>
<div class="col">
<p>This is the combination Proposition 1 appeared to rule out: an influence function that is
unbounded above, so the large observations retain full leverage on the shape parameter, and a
functional that is nonetheless local to the upper tail. It is an inversion of Huber's knot &mdash;
constant influence where Huber is linear, linear influence where Huber is bounded &mdash; and the
motivation is the inversion of Huber's premise.</p>

<p>Measured on the contaminated truth, the claim holds, but only when the comparison is made
carefully. Read at matched <em>nominal</em> level the one-sided functional appears to shed all
contamination at any c, which is too good: at c = 0.1 the nominal p = 0.90 functional sits at the
0.997 quantile, so small c does not shed the contamination but outruns it. Compared at matched
<em>effective</em> level F(t), surviving contamination at the 0.95 effective level is 0.000% for
c &le; 0.5, against 0.569% for the quantile and 15.7% for the expectile. The result survives the
fair comparison.</p>

<div class="note"><p>Two limits are worth recording. As c &rarr; &infin; the expectile is
recovered exactly. As c &rarr; 0 the functional diverges to the upper endpoint, so the family does
not interpolate to the quantile; c is a one-sided dial away from the expectile rather than a bridge
between two usable estimators. A symmetric inversion, &psi;(u) = sign(u) max(|u|, c), does run
expectile to quantile, but sheds contamination only as it converges to the quantile, and is not
pursued further.</p></div>
</div>
"""

S5 = """
<h2><span class="n">5</span><span>Simulation study</span></h2>
<h3><span class="n">5.1</span>Design</h3>
<div class="col">
<p>The truth is the distribution of max(X, Z) with X ~ GEV(0, 1, 0.2) and Z ~ Normal(1.5, 0.8)
independent, so its distribution function is the product F<sub>GEV</sub>&#8202;&Phi;. The normal
dominates the body and vanishes in the tail: the GEV is the <em>exactly correct tail model</em> and
a <em>wrong body model</em>, which is the situation the method is built for. On a return-period
axis the two curves are indistinguishable beyond about T = 15; the damage is below the 90th
percentile, where the truth's median is 4.6 times its GEV component's and its mean 2.4 times.
Asymptotically, maximum likelihood underestimates the 1000-year return level by 39% under this
contamination and L-moments by 32%.</p>

<p>Throughout, n = 100 &mdash; a long annual-maximum record &mdash; with 2000 Monte Carlo
replicates, return periods on a log grid from 2 to 1000 years, and paired standard errors for all
comparisons between estimators, since every estimator sees the same datasets. Where a graft onto an
empirical body is used, the construction is the smooth graft of [<span
style="color:var(--accent-2)">smooth-graft reference</span>], with the handover weight taken equal
to the fitting weight.</p>
</div>

<h3><span class="n">5.2</span>Under misspecification</h3>
<div class="col">
<p>Both composite estimators remove essentially all of the asymptotic bias, cutting 39% and 32% to
under 1%. The question is variance, and this is where the choice of loss matters: <strong>at matched
weight, the composite expectile estimator's mean squared error for the 1000-year return level is
0.17 to 0.58 times the composite quantile estimator's</strong> &mdash; a two- to six-fold
improvement, at every weight tested, with the L1 version's far-tail standard deviation roughly
double the L2 version's. The variance penalty that makes the L1 construction unattractive is
substantially repaired by the change of loss.</p>

<p>Mixing helps further. Sweeping &alpha; in (3) inside a GPD study with a convex weight, the MSE
relative to peaks-over-threshold is 0.74 at T = 50 and 0.79 at T = 100 for &alpha; = 0.5, against
0.89 and 1.15 for pure L1 and 0.86 and 0.87 for pure L2: a gain of up to 14% over the better of the
two pure losses, tapering to nothing beyond T = 500 where L2 already wins. The shape parameter is
also best recovered near &alpha; = 0.5. The best &alpha; rises with the return period, and an oracle
that chooses it separately at every level improves on the best single &alpha; by 14% at the far tail
and 4.4% on average &mdash; the entire prize available to a level-dependent schedule, before any
cost of estimating one.</p>

<p>The one-sided influence function of Proposition 2 does better still. With the knot set as
c = k&#8202;&times;&#8202;IQR(y), fixed from the data before optimising so that &rho; does not move
with &theta;, k = 4 beats pure L2 at every return period from 107 to 1000, with paired t statistics
between &minus;5.8 and &minus;12.5, and reaches 0.179 relative to the reference at T = 1000 against
a previous best of 0.209. The optimum in k is interior &mdash; k &rarr; &infin; is pure L2 exactly,
and the curve turns back up on both sides of k = 4 (Figure 2) &mdash; and it is wide rather than a
knife edge, with k = 16 still beating L2 at t = &minus;7.9.</p>

<div class="note"><p>The mechanism is not what produces the gain. Proposition 2 predicts less bias;
finite-sample the bias is <em>worse</em>, &minus;2.55 at T = 1000 for k = 4 against &minus;1.94 for
pure L2. The estimator is more biased and wins on variance. And the regime where the theory helps
most &mdash; c &le; 0.5, where contamination vanishes &mdash; is catastrophic at n = 100, giving MSE
ratios above 4 at short return periods, because 10.1% of the fitting weight then lands above the
largest observation (Section 3.3). The result should be presented as an empirical finding with a
suggestive derivation attached, not as theory confirmed by simulation.</p></div>
</div>
FIG_KNOT

<h3><span class="n">5.3</span>Under correct specification</h3>
<div class="col">
<p>None of the above is deployable without knowing what it costs when the family is right. We
therefore repeat the study on data drawn from the fitted family exactly. On iid GEV(0, 1, 0.2) at
n = 100, where maximum likelihood is efficient, the composite estimators lose (Table 2).</p>
</div>
<div class="tablewrap">
<table>
<caption><b>Table 2.</b> MSE relative to the GEV MLE on data that really are GEV, n = 100, 2000
replicates. Above 1 is the price paid. The final column is the median fitted shape against a true
&xi; = 0.2.</caption>
<thead><tr><th>estimator</th><th>T=2</th><th>T=10</th><th>T=48</th><th>T=203</th><th>T=529</th><th>T=1000</th><th>median &xi;&#770;</th></tr></thead>
<tbody>
<tr class="ref"><td>GEV L-moments</td><td>1.05</td><td>0.95</td><td>0.94</td><td>0.98</td><td>1.01</td><td>1.03</td><td>0.184</td></tr>
<tr><td>composite L1</td><td class="lose">3.68</td><td>1.20</td><td>1.11</td><td class="lose">1.64</td><td class="lose">2.24</td><td class="lose">2.83</td><td>0.167</td></tr>
<tr><td>composite L2</td><td class="lose">2.03</td><td class="lose">2.03</td><td>1.30</td><td>1.28</td><td>1.31</td><td>1.32</td><td class="lose">0.095</td></tr>
<tr><td>&alpha;-elastile, &alpha; = 0.5</td><td class="lose">2.01</td><td>1.51</td><td class="win">0.98</td><td>1.04</td><td>1.15</td><td>1.25</td><td>0.122</td></tr>
<tr><td>one-sided, k = 4</td><td class="lose">2.47</td><td class="lose">2.28</td><td>1.43</td><td>1.26</td><td>1.18</td><td class="win">1.12</td><td class="lose">0.085</td></tr>
<tr class="ref"><td>GEV MLE</td><td>1.00</td><td>1.00</td><td>1.00</td><td>1.00</td><td>1.00</td><td>1.00</td><td>0.198</td></tr>
</tbody></table>
</div>
<div class="col">
<p>Three readings. The body is where the loss hurts: every composite estimator is two to nearly four
times the MLE's error at T = 2, the weight doing precisely what it was told in a setting where the
body was informative. The far-tail premium is modest for two members &mdash; 12% for the one-sided
loss and 25% for the elastile &mdash; which is a defensible insurance cost. And the composite
quantile estimator gets <em>worse</em> with return period rather than better, 1.11 at T = 48 rising
to 2.83 at T = 1000: its variance penalty appears with no misspecification to justify it, so the L1
construction has no regime in which it is preferred.</p>

<div class="note"><p><b>The bias runs in the unsafe direction.</b> At T = 1000 the MLE is biased
+1.01 and L-moments +0.45, while the composite expectile is &minus;1.52 and the one-sided loss
&minus;1.93. The composite estimators under-predict the far tail, which for design purposes is the
dangerous direction. The cause is in the last column of Table 2: tail-weighted fitting underestimates
the shape parameter badly at n = 100 even with a perfect family, and a low shape shrinks the
extrapolation.</p></div>

<p>The GPD case adds a caution about how return periods are read. Fitting a GPD with the threshold
known at zero &mdash; two parameters, the same two for every estimator, a regular problem in which
the MLE is efficient &mdash; L-moments beats the MLE at every return period, the known small-sample
result, so it is the honest yardstick. Against it, only the one-sided loss wins, and only beyond
about T = 26 (Table 3). Crucially, when a GPD is deployed above a threshold with exceedance rate
&zeta;, the T-year level solves &zeta;&#8202;S(x) = 1/T, so a column labelled T here is the
(T/&zeta;)-year level. With &zeta; = 0.1 the 50-to-1000-year range that design work occupies is
T = 5 to 100 in this table &mdash; the middle, not the right edge.</p>
</div>
<div class="tablewrap">
<table>
<caption><b>Table 3.</b> MSE relative to GPD L-moments on data that really are GPD, threshold known,
n = 100. All rows grafted. Below 1 is better than L-moments. The second row gives the corresponding
return period when the GPD sits above a threshold with exceedance rate &zeta; = 0.1.</caption>
<thead><tr><th>estimator</th><th>T=2</th><th>T=5</th><th>T=10</th><th>T=26</th><th>T=48</th><th>T=107</th><th>T=1000</th></tr></thead>
<tbody>
<tr class="sub"><td>&zeta; = 0.1: true return period</td><td>20</td><td>52</td><td>98</td><td>256</td><td>484</td><td>1074</td><td>10<sup>4</sup></td></tr>
<tr class="ref"><td>GPD MLE</td><td>1.06</td><td>1.00</td><td>1.03</td><td>1.09</td><td>1.12</td><td>1.15</td><td>1.20</td></tr>
<tr class="ref"><td>empirical distribution alone</td><td>1.43</td><td>1.51</td><td>1.59</td><td>1.55</td><td>1.44</td><td>1.33</td><td>1.05</td></tr>
<tr><td>composite L1</td><td class="lose">1.36</td><td class="lose">1.36</td><td class="lose">1.26</td><td class="lose">1.21</td><td class="lose">1.25</td><td class="lose">1.32</td><td class="lose">1.51</td></tr>
<tr><td>composite L2</td><td class="lose">1.36</td><td class="lose">1.35</td><td class="lose">1.19</td><td>1.00</td><td>0.98</td><td>1.00</td><td>0.97</td></tr>
<tr><td>&alpha;-elastile, &alpha; = 0.5</td><td class="lose">1.35</td><td class="lose">1.36</td><td class="lose">1.21</td><td>1.03</td><td>1.02</td><td>1.03</td><td>1.01</td></tr>
<tr><td>one-sided, k = 4</td><td class="lose">1.36</td><td class="lose">1.35</td><td class="lose">1.19</td><td>0.98</td><td class="win">0.94</td><td class="win">0.91</td><td class="win">0.76</td></tr>
</tbody></table>
</div>
<div class="col">
<p>Across the practically occupied range, then, the best composite estimator is 35% worse than
L-moments at the 50-year level, 19% worse at 100, level at 250, and 9% better at 1000. The far-tail
advantage is real but sits largely beyond the horizon of interest.</p>
</div>
FIG_CORRECT

<h3><span class="n">5.4</span>Decomposing the apparent gain</h3>
<div class="col">
<p>The correctly-specified GPD study also exposes a confound that we believe is general. Compared
against peaks-over-threshold with maximum likelihood at a 0.90 threshold &mdash; the conventional
reference &mdash; the composite estimators win comfortably in the far tail on the contaminated
truth, at an MSE ratio of 0.209 at T = 1000. But they also win with <em>no misspecification at
all</em>, at a ratio of 0.344. The reason is that at n = 100 a 0.90 threshold leaves ten exceedances
to fit two parameters, while the composite estimators use all one hundred observations; measured
against an efficient full-sample estimator, POT&ndash;MLE(0.90) is itself 2.58 times worse at
T = 1000 despite being correctly specified.</p>
</div>
<div class="prop"><span class="lab">The decomposition</span>
<p>On a log scale, <strong>68% of the apparent gain over the conventional reference is present with
a perfectly specified model</strong>, and only 32% is attributable to the misspecification the
method exists to exploit. The gain is real, but it is mostly data efficiency, and a comparison
reported only on a misspecified truth does not distinguish the two.</p></div>
<div class="col">
<p>The remedy is cheap: report the same comparison on data simulated from the fitted family, and
include a reference that does not discard observations. We would suggest this as routine for any
proposal that claims to improve on peaks-over-threshold, since the threshold is doing work that is
easily mistaken for the method's.</p>

<p>A related point concerns the shortest return periods. When the composite fit is grafted onto an
empirical body, the T = 2 error is 1.284, 1.288, 1.283 and 1.289 for four very different tail
models, and the empirical distribution alone gives 1.353. Four indistinguishable numbers indicate
that the body is doing all the work: that column measures the empirical distribution against a
parametric fit, and is a cost of grafting rather than of the loss.</p>
</div>
"""

S6 = """
<h2><span class="n">6</span><span>Discussion</span></h2>
<div class="col">
<p>Tail-weighted composite losses do what they are designed to do. Under body misspecification of
the kind that rain-on-snow catchments produce, they remove essentially all of the asymptotic bias
that likelihood and L-moments incur, and the change from the pinball loss to the asymmetric square
removes most of the variance penalty that made the L1 construction unattractive. The one-sided
influence function of Proposition 2 is the most accurate far-tail member we examined, and it has
the distinction of being derived from the obstacle rather than interpolated between two existing
losses.</p>

<p>Three qualifications should travel with that conclusion.</p>

<p>First, the method is insurance and should be priced as such. On a correctly specified GEV the
composite estimators cost 12&ndash;32% in the far tail and two to four times in the body; on a
correctly specified GPD they cost 20&ndash;35% across the return periods design work occupies, and
repay only beyond roughly the 250-year level. The question a practitioner must answer first is
therefore not which loss to use but whether the body model is wrong &mdash; and the argument that
it usually is, for a catchment with two flood-generating mechanisms, is a physical argument rather
than a statistical one.</p>

<p>Second, the errors run the wrong way for design. Every composite estimator we examined
under-predicts the far tail when the model is correct, because tail-weighted fitting underestimates
the shape parameter at realistic sample sizes. Where the consequence of under-prediction is
asymmetric &mdash; as it is for flood defence &mdash; this may outweigh a modest MSE advantage, and
an explicit bias correction or a shrinkage-aware reporting convention would be worth developing.</p>

<p>Third, and most consequential for the literature, the conventional comparison overstates the
case. Two thirds of the apparent advantage over peaks-over-threshold in our misspecified study is
present when nothing is misspecified, because the reference discards nine tenths of the sample. We
do not think this is peculiar to our estimators; any whole-sample method compared against a
high-threshold POT fit will inherit the same advantage, and the practice of reporting only the
misspecified comparison makes it invisible.</p>

<p>Several directions remain. The exponent of the loss can be varied continuously between one and
two, giving the L<sup>a</sup>-quantiles; we found these dominated by the elastile on accuracy at
every return period tested, but they are the only continuous path off the pinball loss that widens
the admissible domain (Table 1), and are the natural recourse when &xi; is too large for the
asymmetric square to have a limit distribution. The mixing weight &alpha; may be made a function of
the level while remaining proper, though the oracle ceiling for doing so is small. And the sandwich
variance of Section 2.2 gives analytic standard errors that we have not exploited; the near-flat
ridge that dominates the estimators' variance is the information matrix approaching singularity
over the tail-weighted region, which suggests a reparameterisation or a penalty in that direction
rather than a further change of loss.</p>
</div>
"""

REFS = """
<h2><span class="n">7</span><span>References</span></h2>
<ol class="refs">
<li><b>Basu, A., Harris, I. R., Hjort, N. L. and Jones, M. C.</b> (1998). Robust and efficient
estimation by minimising a density power divergence. <i>Biometrika</i> <b>85</b>, 549&ndash;559.</li>
<li><b>Bradic, J., Fan, J. and Wang, W.</b> (2011). Penalized composite quasi-likelihood for
ultrahigh dimensional variable selection. <i>Journal of the Royal Statistical Society, Series B</i>
<b>73</b>, 325&ndash;349.</li>
<li><b>Breckling, J. and Chambers, R.</b> (1988). M-quantiles. <i>Biometrika</i> <b>75</b>,
761&ndash;771.</li>
<li><b>Dupuis, D. J.</b> (1998). Exceedances over high thresholds: a guide to threshold selection.
<i>Extremes</i> <b>1</b>, 251&ndash;261.</li>
<li><b>Fung, T. C.</b> (2022). Maximum weighted likelihood estimator for robust heavy-tail
modelling of finite mixture models. <i>Insurance: Mathematics and Economics</i> <b>107</b>,
180&ndash;198.</li>
<li><b>Gomes, M. I., de Haan, L. and Rodrigues, L. H.</b> (2008). Tail index estimation for
heavy-tailed models: accommodation of bias in weighted log-excesses. <i>Journal of the Royal
Statistical Society, Series B</i> <b>70</b>, 31&ndash;52.</li>
<li><b>Hosking, J. R. M. and Wallis, J. R.</b> (1987). Parameter and quantile estimation for the
generalized Pareto distribution. <i>Technometrics</i> <b>29</b>, 339&ndash;349.</li>
<li><b>Huber, P. J.</b> (1964). Robust estimation of a location parameter. <i>Annals of
Mathematical Statistics</i> <b>35</b>, 73&ndash;101.</li>
<li><b>Ju&aacute;rez, S. F. and Schucany, W. R.</b> (2004). Robust and efficient estimation for the
generalized Pareto distribution. <i>Extremes</i> <b>7</b>, 237&ndash;251.</li>
<li><b>Koenker, R. and Bassett, G.</b> (1978). Regression quantiles. <i>Econometrica</i> <b>46</b>,
33&ndash;50.</li>
<li><b>Newey, W. K. and Powell, J. L.</b> (1987). Asymmetric least squares estimation and testing.
<i>Econometrica</i> <b>55</b>, 819&ndash;847.</li>
<li><b>Zou, H. and Yuan, M.</b> (2008). Composite quantile regression and the oracle model
selection theory. <i>Annals of Statistics</i> <b>36</b>, 1108&ndash;1126.</li>
</ol>

<footer>
<p><b>Draft note, not part of the manuscript.</b> Bibliographic details above are from memory and a
literature search and should be verified against the sources before submission; the smooth-graft
citation in Section 5.1 is a placeholder. All numerical results are reproducible from the
accompanying simulation code. Suggested outlets, in order: <i>Extremes</i> (the theory in Sections
3 and 4 carries the paper), <i>Environmetrics</i>, or <i>Water Resources Research</i> if the
rain-on-snow application is developed into a case study.</p>
</footer>
</main>
</div>
"""


def build():
    return "".join([
        HEAD, S1,
        S2.replace("FIG_PSI", fig(
            "fig-loss-geometry.png", 1,
            "<b>The family, and why robust statistics passed it over.</b> Left, loss balls in "
            "residual space: the set of residual vectors with &Sigma;&rho;(r<sub>i</sub>) &le; 1 "
            "is a cross-polytope for the pinball loss and a ball for the asymmetric square, with "
            "the elastile between. Right, influence functions. The elastile's &psi; is unbounded "
            "for every &alpha; &gt; 0, so its gross-error sensitivity is infinite and its "
            "breakdown point zero; Huber's is bounded. The one-sided inversion of Proposition 2 "
            "follows Huber below the knot and the expectile above it, which is the combination "
            "the tail problem calls for.")),
        S3,
        S4,
        S5.replace("FIG_KNOT", fig(
            "fig-invhuber-gpd.png", 2,
            "<b>The knot sweep under misspecification.</b> Left, MSE relative to "
            "peaks-over-threshold against return period for knots k = 0.25 to 16, with pure L2 "
            "(the k &rarr; &infin; limit) dashed and the &alpha; = 0.5 elastile dotted. Small "
            "knots are catastrophic below T = 200, the mechanism being Section 3.3. Right, the "
            "same numbers read along k: the optimum is interior at k = 4, circled, and the curve "
            "turns back towards the pure-L2 line on both sides.")).replace("FIG_CORRECT", fig(
            "fig-correct-gpd2.png", 3,
            "<b>The premium, on a correctly specified GPD.</b> Threshold known, so every method "
            "estimates the same two parameters. Left, ungrafted: the composite fits are 1.5 to "
            "2.0 times the MLE at T = 2 because the weight neglected the body. Right, grafted "
            "onto an empirical body, which brings all four to about 1.29 there &mdash; the same "
            "value for four different tail models, close to the empirical distribution's own 1.35 "
            "(dotted), because the body is doing the work. L-moments, dashed, beats the MLE "
            "throughout.")),
        S6, REFS,
    ])


if __name__ == "__main__":
    html = build()
    (ROOT / "paper" / "manuscript.html").write_text(html)
    print("wrote paper/manuscript.html  (%.2f MB)" % (len(html) / 1e6))
