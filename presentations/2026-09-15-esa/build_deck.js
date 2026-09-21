// Build the ESA progress-meeting deck (15 September 2026) as an editable .pptx.
//
//   npm install pptxgenjs@3.12.0      (anywhere on NODE_PATH)
//   node presentations/2026-09-15-esa/build_deck.js
//
// pptxgenjs rather than python-pptx because Keynote imports its output (the
// internal 2026-09-10 deck was made with it) and refused the python-pptx file.
// Styling follows that internal deck: Cambria titles, Calibri body, navy section
// slides, orange kickers. Figures in figs/ come from it. The Gantt chart is a
// native table transcribed from `new gantt chart.docx`, so it can be edited.
const path = require("path");
const pptxgen = require("pptxgenjs");

const HERE = __dirname;
const FIG = (f) => path.join(HERE, "figs", f);
const OUT = path.join(HERE, "2026-09-15 ESA progress meeting.pptx");

const NAVY = "16333F", ORANGE = "C1440E", LABEL = "1F5673", GREY = "5A6A72";
const PINK = "F6E3DC", MIST = "EDF2F5", STEEL = "8FB3C4", PALE = "9DBBCA";
const RULE = "D5DEE3", WHITE = "FFFFFF", FOCUS = "1F5673", ONGOING = "C6D8E1";
const SERIF = "Cambria", SANS = "Calibri";

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10 x 5.625 in
pres.author = "Vincenzo Coia";
pres.title = "Rain-on-snow: ESA progress meeting";

// ---- helpers ---------------------------------------------------------------

// runs: a string, or an array of paragraphs; each paragraph is a string or an
// array of [text, overrides] pairs.
function text(slide, x, y, w, h, runs, o = {}) {
  const base = {
    fontFace: o.font || SANS, fontSize: o.size || 15, color: o.color || NAVY,
    bold: !!o.bold, italic: !!o.italic,
  };
  const paras = typeof runs === "string" ? [runs] : runs;
  const items = [];
  paras.forEach((para, i) => {
    const pieces = typeof para === "string" ? [[para, {}]] : para;
    pieces.forEach(([t, over], j) => {
      const opts = {
        fontFace: over.font || base.fontFace, fontSize: over.size || base.fontSize,
        color: over.color || base.color,
        bold: over.bold !== undefined ? over.bold : base.bold,
        italic: over.italic !== undefined ? over.italic : base.italic,
      };
      if (j === pieces.length - 1 && i < paras.length - 1) {
        opts.breakLine = true;
        if (o.spaceAfter) opts.paraSpaceAfter = o.spaceAfter;
      }
      items.push({ text: t, options: opts });
    });
  });
  slide.addText(items, {
    x, y, w, h, margin: 0, valign: o.valign || "top", align: o.align || "left",
    isTextBox: true, paraSpaceAfter: o.spaceAfter || 0,
  });
}

const rect = (slide, x, y, w, h, fill) =>
  slide.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: fill, width: 0 } });

const hline = (slide, x, y, w, color = RULE, width = 0.75) =>
  slide.addShape(pres.shapes.LINE, { x, y, w, h: 0, line: { color, width } });

const vline = (slide, x, y, h, color, width) =>
  slide.addShape(pres.shapes.LINE, { x, y, w: 0, h, line: { color, width } });

function header(slide, kicker, title, subtitle) {
  text(slide, 0.5, 0.28, 9, 0.24, kicker.toUpperCase(), { size: 11, color: ORANGE, bold: true });
  text(slide, 0.5, 0.52, 9, 0.56, title, { size: 28, font: SERIF, bold: true });
  if (subtitle) text(slide, 0.5, 1.11, 9, 0.30, subtitle, { size: 13.5, color: GREY });
}

function section(word, num, title, blurb) {
  const s = pres.addSlide();
  s.background = { color: NAVY };
  text(s, 0.5, 1.70, 9, 0.28, `SECTION ${word}`, { size: 12, color: ORANGE, bold: true });
  text(s, 0.5, 2.14, 9, 0.90, title, { size: 38, font: SERIF, bold: true, color: WHITE });
  text(s, 0.5, 3.14, 9, 0.70, blurb, { size: 15, color: PALE });
  text(s, 8.9, 4.66, 0.6, 0.58, String(num), { size: 30, font: SERIF, bold: true, color: "2C4C5C", align: "right" });
  return s;
}

function kvRows(slide, rows, top, rowH, labelW = 2.35, gap = 0.15, size = 15) {
  rows.forEach(([label, value], i) => {
    const y = top + i * rowH;
    if (i > 0) hline(slide, 0.5, y - 0.06, 9);
    text(slide, 0.5, y, labelW, rowH - 0.1, label, { size, color: LABEL, bold: true });
    text(slide, 0.5 + labelW + gap, y, 9 - labelW - gap, rowH - 0.1, value, { size });
  });
}

// ---- 1. Title ----------------------------------------------------------------

let s = pres.addSlide();
s.background = { color: NAVY };
text(s, 0.5, 1.75, 9, 0.85, "Rain-on-snow", { size: 46, font: SERIF, bold: true, color: WHITE });
text(s, 0.5, 2.62, 9, 0.75, "Progress update", { size: 38, font: SERIF, color: STEEL });
text(s, 0.5, 3.62, 9, 0.30, "Vincenzo Coia   ·   ESA progress meeting   ·   15 September 2026", { size: 15, color: "B9CDD8" });
s.addNotes("Four parts: software updates, rain-on-snow findings, statistical developments, and where we are on the timeline.");

// ---- 2. Software ---------------------------------------------------------------

s = pres.addSlide();
header(s, "Updates", "probaverse: new releases on CRAN", "distionary and distplyr both have new versions out.");
rect(s, 0.5, 1.62, 9, 0.78, PINK);
text(s, 0.78, 1.62, 8.44, 0.78,
  "The graft machinery is now published in the packages, rather than living as loose code in the project repository.",
  { size: 15, valign: "middle" });
kvRows(s, [
  ["Supports", "A distribution knows where it lives, not only what it evaluates to. Now complete."],
  ["Discrete parts", "Algorithms anticipate distributions with a discrete component, and are faster for it."],
  ["Hard graft", "Supported in distplyr. The smooth graft comes later in this deck."],
  ["JOSS", "Software paper submitted."],
], 2.66, 0.66);
s.addNotes("Why supports matter: grafting a tail onto an empirical body joins something that stops dead to something that runs forever. Without a support, the algorithms cannot tell an impossible value from an unlikely one.");

// ---- 3. Section: rain-on-snow ------------------------------------------------------

section("ONE", 1, "Rain-on-snow", "Which mix of rain and snowmelt sits behind a T-year runoff event.");

// ---- 4. Dials of risk ----------------------------------------------------------

s = pres.addSlide();
header(s, "Rain-on-snow", "Dials of risk", "Rain and snowmelt behind the T-year runoff event, as T grows.  Cell 3.");
s.addImage({ path: FIG("dials-of-risk.gif"), x: 0.45, y: 1.59, w: 4.90, h: 3.675 });
text(s, 5.65, 1.58, 3.85, 0.28, "What to watch", { size: 12, color: LABEL, bold: true });
text(s, 5.65, 1.94, 3.85, 0.80,
  "As T grows, the likeliest mix moves towards rain alone: more rainfall, with an ordinary snowmelt.", { size: 14 });
rect(s, 5.65, 2.86, 3.85, 2.34, PINK);
text(s, 5.90, 2.98, 3.36, 2.14, [
  [["Not yet reliable.", { bold: true, color: ORANGE, size: 15 }]],
  "Key problems in the distributional learning model sit under this picture: the learned conditionals stop where the data stop, so the rare-event end is the least trustworthy part.",
], { size: 13.5, spaceAfter: 6 });
s.addNotes("At T = 2 the likeliest mix is a 2.1-year rainfall with a 0.3-year snowmelt; at T = 200, a 9.4-year rainfall with a 0.24-year snowmelt. Rainfall explains over 90% of the variation in runoff at peaks in this cell. The animation plays in slideshow mode.");

// ---- 5. The next question ------------------------------------------------------------

s = pres.addSlide();
header(s, "Rain-on-snow", "The next question", "Condition on the state of the snowpack, not only on what is falling on it.");
rect(s, 0.5, 1.62, 9, 1.22, PINK);
text(s, 0.85, 1.62, 8.3, 1.22,
  "For a T-year runoff event, with a given amount of water available in the snowpack, how much rain is enough to trigger it?",
  { size: 21, font: SERIF, bold: true, valign: "middle" });
kvRows(s, [
  ["The new driver", "Water available in the snowpack, which is closely tied to snow wetness."],
  ["Where ESA comes in", "Snow wetness is observed from space: the natural moment to bring in an ESA data product."],
  ["What is ready", "The analysis pipeline now takes any driver and answers this question directly."],
], 3.12, 0.72);
s.addNotes("The pipeline was restructured so each analysis is a choice of drivers, a learning model, and a question. A first pass with snowmelt standing in for available water already runs; the forest used there cannot yet resolve events rarer than about 2 years, which is the same tail problem as the previous slide.");

// ---- 6. Section: statistics ------------------------------------------------------

section("TWO", 2, "Statistical developments", "Tail estimation now pays off. The learning model has two new challenges.");

// ---- 7. The smooth graft, revisited --------------------------------------------------

s = pres.addSlide();
header(s, "Tail estimation", "The smooth graft, revisited");
[
  [0.5, MIST, LABEL, "Where we left off", "The smooth graft on its own did not look like it was worth anything.", "A gentler handover of the same tail."],
  [5.15, PINK, ORANGE, "Where we are now", "Combined with a new composite estimation method, built on expectiles, it improves on classical tail estimation.", "The composite weight doubles as the handover."],
].forEach(([x, fill, tc, title, body, foot]) => {
  rect(s, x, 1.45, 4.35, 2.45, fill);
  text(s, x + 0.3, 1.67, 3.78, 2.05, [
    [[title, { bold: true, color: tc, size: 17 }]],
    body,
    [[foot, { italic: true, color: GREY, size: 13 }]],
  ], { size: 15, spaceAfter: 10 });
});
rect(s, 0.5, 4.15, 9, 0.85, MIST);
text(s, 0.78, 4.15, 8.44, 0.85, [
  [["The composite estimator: ", { bold: true, color: LABEL }], ["an expectile loss integrated over levels, weighted towards the tail", {}]],
  [["S(θ)  =  Σᵢ  ∫₀¹  w(p) · | p − 1(yᵢ < E(p|θ)) | · (yᵢ − E(p|θ))²  dp", { font: SERIF, size: 13 }]],
], { size: 12.5, valign: "middle", spaceAfter: 3 });
s.addNotes("E is the family's expectile function. The weight w(p) does the work; the useful range puts it low, p0 around 0.5 to 0.8, not 0.95.");

// ---- 8. Against practitioners ------------------------------------------------------------

s = pres.addSlide();
header(s, "Tail estimation", "Against what practitioners do", "Composite fits: smooth grafted.  Reference: POT-MLE at the 90th percentile, hard grafted.");
s.addImage({ path: FIG("tail-vs-practitioners.png"), x: 0.6, y: 1.46, w: 8.8, h: 3.52 });
text(s, 0.5, 5.04, 9, 0.50,
  "Composite L2, w = p⁶, smooth grafted: MSE ratio 0.86 at T = 50 falling to 0.21 at T = 1000; closer on 77% of datasets.",
  { size: 14, italic: true, color: ORANGE });
s.addNotes("Threshold choice dominates the competitor: POT-MLE at 0.85 / 0.90 / 0.95 gives MSE in the ratio 1.62 : 1.00 : 0.52 at T = 1000, a threefold spread from a choice the data cannot make. The composite estimator has no threshold; the weight plays that role, smoothly.");

// ---- 9. Challenge 1 -----------------------------------------------------------------

s = pres.addSlide();
header(s, "Distributional learning", "Challenge 1: the marginal of runoff", "Averaging the learned conditionals over the drivers is ill-conditioned, for two reasons.");
s.addImage({ path: FIG("marginal-by-mixing.png"), x: 0.45, y: 1.62, w: 5.55, h: 5.55 * 780 / 1700 });
text(s, 6.2, 1.62, 3.3, 3.0, [
  [["(i)  The maximum of noise", { bold: true, color: LABEL, size: 15 }]],
  "The mixture's tail index is the largest tail index among the conditionals: the maximum of a set of noisy estimates.",
  [["(ii)  The jump needs extrapolation", { bold: true, color: LABEL, size: 15 }]],
  "The marginal EVI is often strictly larger than every conditional EVI. Reaching it takes predictions far out in predictor space, and machine learning cannot go there.",
], { size: 13.5, spaceAfter: 6 });
rect(s, 0.5, 4.55, 9, 0.62, PINK);
text(s, 0.78, 4.55, 8.44, 0.62, "Learning stops at the data.", { size: 17, font: SERIF, bold: true, color: ORANGE, valign: "middle" });
s.addNotes("Point (ii) is from Coia et al. (2020): the marginal EVI is always at least the largest conditional EVI, strictly so for many copulas. Closing that gap needs the conditionals to hold indefinitely into the tail of x.");

// ---- 10. Challenge 2 ----------------------------------------------------------------

s = pres.addSlide();
header(s, "Distributional learning", "Challenge 2: copula transport", "Like peaks-over-threshold: the copula carries learned distributions beyond the predictor space.");
s.addImage({ path: FIG("copula-transport.png"), x: 1.09, y: 1.48, w: 7.83, h: 7.83 * 760 / 1750 });
text(s, 0.5, 4.95, 9, 0.55, [
  [["Left: it reproduces a conditional we can check. Right: transports pushed past the data.  ", {}],
   ["Still open: when does it give favourable estimates?", { bold: true }]],
], { size: 14, italic: true, color: ORANGE });
s.addNotes("One conditional plus the copula determines the whole family. Combine several by the median survival, not the mean, which would propagate the most extreme tail index again. So far the noise in distributional learning has struggled to beat estimating the distribution of runoff directly; that is what needs testing.");

// ---- 11. Papers -------------------------------------------------------------------

s = pres.addSlide();
header(s, "Statistical developments", "Three papers", "One for each piece of the method, at different stages.");
[
  ["Smooth graft", "When and how to hand off to a tail model.",
   "Loose ends in the simulation section, then ready for final review.", "Nearly ready", NAVY],
  ["Composite estimation", "Improved tail estimation.",
   "The bones are in place, and a compelling simulation study is already there.", "Drafting", LABEL],
  ["Copula transport", "Extending learned distributions outside the hull of x.",
   "Inception. The need for it has been clearly demonstrated.", "Inception", ORANGE],
].forEach(([title, purpose, status, chip, chipColor], i) => {
  const x = 0.5 + i * 3.08;
  rect(s, x, 1.62, 2.84, 3.72, i === 2 ? PINK : MIST);
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, {
    x: x + 0.25, y: 1.85, w: 1.25, h: 0.28, rectRadius: 0.14,
    fill: { color: chipColor }, line: { color: chipColor, width: 0 },
  });
  text(s, x + 0.25, 1.85, 1.25, 0.28, chip, { size: 10.5, bold: true, color: WHITE, align: "center", valign: "middle" });
  text(s, x + 0.25, 2.3, 2.4, 0.5, title, { size: 16, font: SERIF, bold: true });
  text(s, x + 0.25, 2.9, 2.34, 0.26, "Purpose", { size: 11, bold: true, color: LABEL });
  text(s, x + 0.25, 3.16, 2.34, 0.75, purpose, { size: 14 });
  text(s, x + 0.25, 3.98, 2.34, 0.26, "Status", { size: 11, bold: true, color: LABEL });
  text(s, x + 0.25, 4.24, 2.34, 0.95, status, { size: 13.5 });
});
s.addNotes("Smooth graft comes first and stands on its own. Composite estimation builds on it. Copula transport addresses the distributional learning challenges on the previous two slides.");

// ---- 12. Gantt -------------------------------------------------------------------

s = pres.addSlide();
header(s, "Project plan", "One year in, and on track", "Phase 2 iterates between statistical development and data, for process insight.");

const TASKS = [ // group, id, name, cells (D focus, L ongoing, . none)
  ["Overhead", "1.1", "Project management", "LLLLLLLLLLLLLLLLLLLLLLLL"],
  ["Overhead", "1.2", "Community", ".........LLLLLLLLLLLLLLL"],
  ["Data", "2.1", "Data compilation", "DD.........DD....D......"],
  ["Data", "2.2", "Event granularity", "...........DD..........."],
  ["Modelling", "3.1", "Min. viable model", ".DD....DD..............."],
  ["Modelling", "3.2", "Extremal learning", ".........DD..DDDD......."],
  ["Modelling", "3.3", "Climate change", ".................DDDDDD."],
  ["Outputs", "4.1", "Software", "..DDDDD..LLLLLLLLLLLLLL."],
  ["Outputs", "4.2", "Flood map insights", "..D....DDLLLLLLLLDDDDDD."],
  ["Outputs", "4.3", "Writing", "......L..LLLLLLLLDDDDDDD"],
];
const DELIV = { "1.1-1": "D1", "4.1-3": "D6", "4.1-7": "D6", "4.2-9": "D3", "4.3-7": "D5",
  "4.3-15": "D2", "4.3-17": "D4", "4.2-23": "D7", "4.3-24": "D8" };
const PHASES = [["Phase 1: MVP + software", 9], ["Phase 2: Iterative development", 8], ["Phase 3: Climate change", 7]];

const gx = 0.5, gy = 1.55, wGroup = 0.78, wId = 0.30, wTask = 1.30;
const gridX = gx + wGroup + wId + wTask;
const cellW = (9.5 - gridX) / 24;
const hPhase = 0.24, hMonth = 0.20, hRow = 0.232;
const tblH = hPhase + hMonth + hRow * TASKS.length;
const border = { type: "solid", pt: 0.5, color: RULE };
const cell = (t, o = {}) => ({
  text: t,
  options: Object.assign({
    fontFace: SANS, fontSize: 8, color: NAVY, fill: { color: WHITE }, border: [border, border, border, border],
    margin: [0, 2, 0, 2], valign: "middle", align: "left",
  }, o),
});

const rows = [];
// Header row 1: blank corner (spans 2 rows x 3 cols) + phases
const r0 = [cell("", { colspan: 3, rowspan: 2, fill: { color: WHITE } })];
PHASES.forEach(([label, n]) => r0.push(cell(label, { colspan: n, fill: { color: MIST }, bold: true, fontSize: 9, align: "center" })));
rows.push(r0);
const r1 = [];
for (let m = 1; m <= 24; m++) r1.push(cell(String(m), { fill: { color: MIST }, fontSize: 7, color: GREY, align: "center" }));
rows.push(r1);
let prevGroup = null;
TASKS.forEach(([group, id, name, cells]) => {
  const r = [];
  if (group !== prevGroup) {
    const n = TASKS.filter((t) => t[0] === group).length;
    r.push(cell(group, { rowspan: n, fill: { color: MIST }, bold: true, fontSize: 8.5, color: LABEL }));
    prevGroup = group;
  }
  r.push(cell(id, { fill: { color: MIST }, fontSize: 7.5, color: GREY }));
  r.push(cell(name, { fill: { color: MIST } }));
  [...cells].forEach((code, k) => {
    const m = k + 1;
    const fill = code === "D" ? FOCUS : code === "L" ? ONGOING : WHITE;
    const d = DELIV[`${id}-${m}`] || "";
    r.push(cell(d, { fill: { color: fill }, fontSize: 6.5, bold: true, align: "center", color: code === "D" ? WHITE : NAVY }));
  });
  rows.push(r);
});
s.addTable(rows, {
  x: gx, y: gy, w: 9.0,
  colW: [wGroup, wId, wTask].concat(Array(24).fill(cellW)),
  rowH: [hPhase, hMonth].concat(Array(TASKS.length).fill(hRow)),
});

const yTop = gy, yBot = gy + tblH;
[9, 17].forEach((mEnd) => vline(s, gridX + cellW * mEnd, yTop, tblH, NAVY, 1.5));
const xNow = gridX + cellW * 12;
vline(s, xNow, yTop + hPhase, tblH - hPhase + 0.12, ORANGE, 2.5);
text(s, xNow - 0.9, yBot + 0.13, 1.8, 0.22, "Today: month 12", { size: 10, bold: true, color: ORANGE, align: "center" });

const ly = yBot + 0.16;
rect(s, 0.5, ly + 0.04, 0.22, 0.13, FOCUS);
text(s, 0.78, ly, 1.1, 0.22, "Focus", { size: 10, color: GREY });
rect(s, 1.45, ly + 0.04, 0.22, 0.13, ONGOING);
text(s, 1.73, ly, 1.3, 0.22, "Ongoing", { size: 10, color: GREY });
text(s, 2.55, ly, 2.2, 0.22, "D1–D8: deliverables", { size: 10, color: GREY });

rect(s, 0.5, 4.72, 9, 0.66, PINK);
text(s, 0.78, 4.72, 8.44, 0.66, [
  [["Next up, as planned: ", { bold: true, color: ORANGE }],
   ["data (2.1–2.2, months 12–13) brings in the snowpack driver; extremal learning (3.2, months 14–17) takes on copula transport.", {}]],
], { size: 13, valign: "middle" });
s.addNotes("Month 1 is October 2025, so month 12 is September 2026. The point: phase 2 alternates statistical development with data and process insight, and the work is where the plan says it should be. The 'possible early end' line from the January replan (after month 7) is no longer shown.");

pres.writeFile({ fileName: OUT }).then((f) => console.log("Wrote " + f));
