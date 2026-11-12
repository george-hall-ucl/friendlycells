// Animated figures for the Algorithm page. Their numbers come from the
// example's results, which friendlycells computes when the site is rendered.
(() => {
  const source = document.getElementById("algorithm-data");
  if (!source) return;
  const data = JSON.parse(source.textContent);
  const NS = "http://www.w3.org/2000/svg";
  // Okabe-Ito, as in plot_neighborhoods()
  const COLORS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00"];
  const labels = data.labels;
  const n = data.samples.length;
  const colorOf = (label) => COLORS[labels.indexOf(label) % COLORS.length];
  const pairs = labels.flatMap((a) => labels.map((b) => [a, b]));
  const short = (name) => name.replace(/^Sample /, "");
  const still = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  const el = (parent, tag, attrs = {}) => {
    const node = document.createElementNS(NS, tag);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
    parent.append(node);
    return node;
  };
  const text = (parent, x, y, content, attrs = {}) => {
    const node = el(parent, "text", {
      x, y, "text-anchor": "middle", "dominant-baseline": "central", ...attrs,
    });
    node.textContent = content;
    return node;
  };
  // A label pair such as A–B, each letter in its label's color
  const pairText = (parent, x, y, [a, b]) => {
    const node = text(parent, x, y, "", { class: "algo-small" });
    for (const [part, color] of [[a, colorOf(a)], ["–", null], [b, colorOf(b)]]) {
      const span = el(node, "tspan", color ? { fill: color, "font-weight": "bold" } : {});
      span.textContent = part;
    }
    return node;
  };
  const sum = (values) => values.reduce((t, v) => t + v, 0);
  const shuffled = (items) => {
    const out = [...items];
    for (let i = out.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [out[i], out[j]] = [out[j], out[i]];
    }
    return out;
  };

  // A round top for an axis, and a step giving a few ticks below it
  const niceScale = (max) => {
    const rough = max / 3;
    const magnitude = 10 ** Math.floor(Math.log10(rough));
    const step = [1, 2, 5, 10].find((m) => m * magnitude >= rough) * magnitude;
    return { top: Math.ceil(max / step - 1e-9) * step, step };
  };
  // A y axis from base (value 0) up height (value top)
  const yAxis = (parent, x, base, height, { top, step }, format) => {
    const g = el(parent, "g", { class: "algo-small" });
    el(g, "line", { x1: x, x2: x, y1: base, y2: base - height, class: "algo-axis" });
    for (let v = 0; v <= top + step / 2; v += step) {
      const y = base - (v / top) * height;
      el(g, "line", { x1: x - 4, x2: x, y1: y, y2: y, class: "algo-axis" });
      text(g, x - 7, y, format(v), { "text-anchor": "end" });
    }
    return g;
  };

  // One playing of a figure. pause() resolves at once when showing the final
  // frame, and rejects once the figure is replayed, ending the old playing.
  const CANCELLED = Symbol("cancelled");
  class Run {
    constructor(animate) {
      this.animate = animate;
      this.cancelled = false;
    }
    pause(ms) {
      if (this.cancelled) return Promise.reject(CANCELLED);
      if (!this.animate) return Promise.resolve();
      return new Promise((resolve, reject) => {
        setTimeout(() => (this.cancelled ? reject(CANCELLED) : resolve()), ms);
      });
    }
  }

  // The largest scale at which every sample fits in width × height
  const fit = (samples, width, height) => Math.min(...samples.map((s) => Math.min(
    width / (Math.max(...s.x) + 1), height / (Math.max(...s.y) + 1),
  )));
  // A sample's links and spots, scaled and centred in the box [x, y, w, h]
  const drawSample = (parent, sample, [x, y, w, h], scale, withLinks) => {
    const dx = x + (w - Math.max(...sample.x) * scale) / 2;
    const dy = y + (h - Math.max(...sample.y) * scale) / 2;
    const g = el(parent, "g", { transform: `translate(${dx} ${dy})` });
    const at = (i) => [sample.x[i] * scale, sample.y[i] * scale];
    const links = withLinks
      ? sample.pairs.map(([i, j]) => {
          const [x1, y1] = at(i);
          const [x2, y2] = at(j);
          return el(g, "line", { x1, y1, x2, y2, class: "algo-link" });
        })
      : [];
    const spots = sample.label.map((label, i) => {
      const [cx, cy] = at(i);
      return el(g, "circle", {
        cx, cy, r: scale * 0.34, fill: colorOf(label), class: "algo-spot",
      });
    });
    return { links, spots };
  };

  // A grid of cells whose shading (0 to 1) and text can be updated
  const cellGrid = (g, size, cell, cls = "") => {
    const cells = Array.from({ length: size }, (_, i) => Array.from({ length: size }, (_, j) => {
      const rect = el(g, "rect", { x: j * cell, y: i * cell, width: cell, height: cell, class: "algo-cell" });
      el(g, "rect", { x: j * cell, y: i * cell, width: cell, height: cell, class: "algo-frame" });
      const value = text(g, j * cell + cell / 2, i * cell + cell / 2, "", cls ? { class: cls } : {});
      return { rect, value };
    }));
    return (i, j, fill, label) => {
      const { rect, value } = cells[i][j];
      rect.style.fillOpacity = 0.85 * fill;
      value.textContent = label;
      value.classList.toggle("algo-on-dark", fill > 0.55);
    };
  };
  // A label × label table
  const drawTable = (parent, x, y, cell, title) => {
    const g = el(parent, "g", { transform: `translate(${x} ${y})` });
    text(g, (labels.length * cell) / 2, -cell * 0.9, title);
    labels.forEach((label, i) => {
      const r = Math.min(7, cell * 0.18);
      el(g, "circle", { cx: i * cell + cell / 2, cy: -cell * 0.35, r, fill: colorOf(label) });
      el(g, "circle", { cx: -cell * 0.35, cy: i * cell + cell / 2, r, fill: colorOf(label) });
    });
    return cellGrid(g, labels.length, cell);
  };
  // A sample × sample matrix
  const drawMatrix = (parent, x, y, cell) => {
    const g = el(parent, "g", { transform: `translate(${x} ${y})` });
    data.names.forEach((name, i) => {
      text(g, i * cell + cell / 2, -12, short(name), { class: "algo-small" });
      text(g, -12, i * cell + cell / 2, short(name), { class: "algo-small" });
    });
    return cellGrid(g, n, cell, "algo-small");
  };

  // The sample comparisons, in order: 1 vs 2, 1 vs 3, 2 vs 3, ...
  const comparisons = data.names.flatMap((_, i) => data.names.map((_, j) => [i, j]))
    .filter(([i, j]) => i < j);
  const squaredDiffs = ([i, j]) => data.profiles[i].map((p, k) => (p - data.profiles[j][k]) ** 2);
  // Each comparison's color, away from the labels' and the samples' blue
  const PAIR_COLORS = ["#D55E00", "#CC79A7", "#7B5EA7"];
  const pairColor = (r) => PAIR_COLORS[r % PAIR_COLORS.length];
  // Each sample's color in the profile bars and comparisons' names
  const SAMPLE_COLORS = ["#0072B2", "#3A9E8F", "#9A9A2E"];
  const sampleColor = (s) => SAMPLE_COLORS[s % SAMPLE_COLORS.length];
  // "Sample 1 vs Sample 2", each sample in its color
  const versus = (node, [i, j], name = (s) => data.names[s]) => {
    node.replaceChildren();
    for (const [part, s] of [[name(i), i], [" vs ", null], [name(j), j]]) {
      const span = el(node, "tspan", s === null ? {} : { fill: sampleColor(s), "font-weight": "bold" });
      span.textContent = part;
    }
  };

  // One row per comparison, shared by both step 5 figures. A row's distance
  // is its squared differences stacked, stretched so the stack's length is
  // their sum's square root.
  const maxD = Math.max(...comparisons.map(([i, j]) => data.distance[i][j]));
  const distanceRows = {
    left: 70,
    width: 290,
    bar: 22,
    draw(svg, top) {
      return comparisons.map(([i, j], r) => {
        const y = top + r * 34;
        const label = text(svg, this.left - 10, y + this.bar / 2, "", {
          "text-anchor": "end", class: "algo-small",
        });
        versus(label, [i, j], (s) => short(data.names[s]));
        const value = text(svg, this.left + this.width + 10, y + this.bar / 2, "", {
          "text-anchor": "start", class: "algo-small",
        });
        return { y, value, label };
      });
    },
    // Where each squared difference sits in the row, at the distance's scale
    segments([i, j]) {
      const d = data.distance[i][j];
      let x = this.left;
      return squaredDiffs([i, j]).map((value) => {
        const w = (value / d ** 2) * (d / maxD) * this.width;
        x += w;
        return { x: x - w, w };
      });
    },
    // A finished row, in its comparison's color
    stack(svg, row, r) {
      const stack = el(svg, "g");
      this.segments(comparisons[r]).forEach(({ x, w }) => {
        el(stack, "rect", { x, y: row.y, width: w, height: this.bar, class: "algo-diff-fill" })
          .style.fill = pairColor(r);
      });
      return stack;
    },
  };

  const figures = {
    // 1. The samples' labelled spots appear, in a random order
    async spots(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 230");
      const width = 640 / n;
      const scale = fit(data.samples, width - 20, 170);
      const spots = data.samples.flatMap((sample, s) => {
        text(svg, s * width + width / 2, 14, data.names[s]);
        return drawSample(svg, sample, [s * width, 30, width, 170], scale, false).spots;
      });
      labels.forEach((label, i) => {
        const x = 320 + (i - (labels.length - 1) / 2) * 60;
        el(svg, "circle", { cx: x - 10, cy: 218, r: 7, fill: colorOf(label) });
        text(svg, x + 6, 218, label);
      });
      for (const spot of spots) spot.classList.add("algo-hidden");
      // In batches, taking a second in all
      const order = shuffled(spots);
      const batches = 25;
      for (let b = 0; b < batches; b++) {
        const from = Math.floor((b * order.length) / batches);
        const to = Math.floor(((b + 1) * order.length) / batches);
        for (const spot of order.slice(from, to)) spot.classList.remove("algo-hidden");
        await run.pause(1000 / batches);
      }
    },

    // 2. Links grow from one spot to the other
    async links(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 300");
      const sample = data.samples[0];
      const scale = fit([sample], 600, 250);
      text(svg, 320, 14, `${data.names[0]}: ${sample.pairs.length} connections`);
      const { links } = drawSample(svg, sample, [20, 40, 600, 250], scale, true);
      if (!run.animate) return;
      for (const link of links) {
        // Grown from either end, at random
        if (Math.random() < 0.5) {
          const [x1, y1] = [link.getAttribute("x1"), link.getAttribute("y1")];
          link.setAttribute("x1", link.getAttribute("x2"));
          link.setAttribute("y1", link.getAttribute("y2"));
          link.setAttribute("x2", x1);
          link.setAttribute("y2", y1);
        }
        const length = Math.hypot(link.getAttribute("x2") - link.getAttribute("x1"),
          link.getAttribute("y2") - link.getAttribute("y1"));
        link.style.strokeDasharray = length;
        link.style.strokeDashoffset = length;
      }
      await run.pause(100);
      // In a random order, in batches, taking 2 seconds in all
      const order = shuffled(links);
      const batches = 40;
      for (let b = 0; b < batches; b++) {
        const from = Math.floor((b * order.length) / batches);
        const to = Math.floor(((b + 1) * order.length) / batches);
        for (const link of order.slice(from, to)) {
          link.classList.add("algo-growing");
          link.style.strokeDashoffset = 0;
        }
        await run.pause(2000 / batches);
      }
    },

    // 3. Each sample's links are counted in turn, into the table below it
    async counts(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 390");
      const width = 640 / n;
      const scale = fit(data.samples, width - 20, 170);
      const cell = 40;
      const max = Math.max(...data.samples.flatMap((s) => s.counts.flat()));
      const panels = data.samples.map((sample, s) => {
        text(svg, s * width + width / 2, 14, data.names[s]);
        const { links } = drawSample(svg, sample, [s * width, 30, width, 170], scale, true);
        const x = s * width + (width - labels.length * cell) / 2 + cell * 0.2;
        const setCell = drawTable(svg, x, 260, cell, "Connection counts");
        const counts = labels.map(() => labels.map(() => 0));
        const show = () => labels.forEach((_, i) => labels.forEach((_, j) => {
          setCell(i, j, counts[i][j] / max, String(counts[i][j]));
        }));
        for (const link of links) link.classList.add("algo-waiting");
        show();
        return { sample, links, counts, show };
      });
      for (const { sample, links, counts, show } of panels) {
        let previous = null;
        for (const [k, [i, j]] of sample.pairs.entries()) {
          const a = labels.indexOf(sample.label[i]);
          const b = labels.indexOf(sample.label[j]);
          // Both directions: A–B and B–A, or A–A twice
          counts[a][b] += 1;
          counts[b][a] += 1;
          previous?.classList.replace("algo-active", "algo-counted");
          links[k].classList.replace("algo-waiting", "algo-active");
          previous = links[k];
          if (run.animate) show();
          await run.pause(5000 / 1.5 / links.length);
        }
        previous?.classList.replace("algo-active", "algo-counted");
        show();
        await run.pause(400);
      }
    },

    // 4. The tables become profiles of counts; a copy slides right and
    // becomes proportions
    async profile(svg, run) {
      const panel = 100;
      const top = 44;
      svg.setAttribute("viewBox", `0 0 640 ${top + panel * n + 16}`);
      const left = 46;
      const shift = 320;
      const step = 28;
      const bar = 20;
      const height = 62;
      const totals = data.samples.map((s) => sum(s.counts.flat()));
      const kinds = {
        counts: {
          scale: niceScale(Math.max(...data.samples.flatMap((s) => s.counts.flat()))),
          title: "Connection counts", values: (s) => data.samples[s].counts.flat(),
          total: (s) => totals[s], format: String,
        },
        shares: {
          scale: niceScale(Math.max(...data.profiles.flat())),
          title: "Proportions of connections", values: (s) => data.profiles[s],
          total: () => 1, format: (v) => v.toFixed(2),
        },
      };
      // A graph drawn at the left, which can be moved and switch kind
      const graph = () => {
        const g = el(svg, "g", { class: "algo-slide" });
        const right = left + 6 + step * pairs.length;
        const title = text(g, (left + right) / 2, 14, "");
        const rows = data.samples.map((_, s) => {
          const base = top + panel * s + 22 + height;
          const total = text(g, left, top + panel * s + 6, "", { "text-anchor": "start", class: "algo-small" });
          el(g, "line", { x1: left, x2: right, y1: base, y2: base, class: "algo-axis" });
          const bars = pairs.map((_, k) => el(g, "rect", {
            x: left + 6 + k * step, y: base - height, width: bar, height, class: "algo-bar",
          }));
          return { base, total, bars, axis: null };
        });
        pairs.forEach((pair, k) => {
          pairText(g, left + 6 + k * step + bar / 2, top + panel * n + 2, pair).classList.add("algo-tiny");
        });
        const show = (kind, grow = true) => {
          const { scale, values, total, format } = kinds[kind];
          title.textContent = kinds[kind].title;
          rows.forEach((row, s) => {
            row.total.textContent = `${data.names[s]}: total ${total(s)}`;
            row.axis?.remove();
            row.axis = yAxis(g, left, row.base, height, scale, format);
            row.bars.forEach((rect, k) => {
              rect.style.transform = grow ? `scaleY(${values(s)[k] / scale.top})` : "scaleY(0)";
            });
          });
        };
        // Each row's count axis rescaled so its top is the proportions'
        // top times the sample's total: the bars reach their proportions
        const rescale = () => {
          const { scale } = kinds.shares;
          title.textContent = "Counts, each axis scaled to its sample's total";
          rows.forEach((row, s) => {
            row.axis?.remove();
            row.axis = yAxis(g, left, row.base, height, scale, (v) => String(Math.round(v * totals[s])));
            row.bars.forEach((rect, k) => {
              rect.style.transform = `scaleY(${data.profiles[s][k] / scale.top})`;
            });
          });
        };
        return { g, show, rescale };
      };
      const counts = graph();
      const copy = graph();
      counts.show("counts");
      copy.show("shares");
      copy.g.style.transform = `translateX(${shift}px)`;
      if (!run.animate) return;

      // The copy matches the counts graph, unseen, before it appears
      copy.g.classList.add("algo-hidden");
      copy.g.style.transform = "none";
      copy.show("counts");
      counts.show("counts", false);
      await run.pause(300);
      counts.show("counts");
      await run.pause(1800);
      copy.g.classList.remove("algo-hidden");
      await run.pause(400);
      copy.g.style.transform = `translateX(${shift}px)`;
      await run.pause(1600);
      copy.rescale();
      await run.pause(2200);
      copy.show("shares");
    },

    // 5. All profiles side by side; for each pair of samples in turn, their
    // squared differences are stacked, and the stack kept as their distance
    async compare(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 370");
      const left = 56;
      const step = 64;
      const bar = 15;
      const base = 210;
      const height = 150;
      const scale = niceScale(Math.max(...data.profiles.flat()));
      const heading = text(svg, 320, 14, "");
      yAxis(svg, left, base, height, scale, (v) => v.toFixed(2));
      el(svg, "line", { x1: left, x2: left + step * pairs.length, y1: base, y2: base, class: "algo-axis" });
      const keys = data.names.map((name, s) => {
        const x = 110 + s * 110;
        const g = el(svg, "g");
        el(g, "rect", { x, y: 30, width: 12, height: 12, fill: sampleColor(s) });
        text(g, x + 18, 36, name, { "text-anchor": "start", class: "algo-small" });
        return g;
      });
      const key = el(svg, "rect", { x: 110 + n * 110, y: 30, width: 12, height: 12, class: "algo-diff-fill" });
      text(svg, 128 + n * 110, 36, "difference", { "text-anchor": "start", class: "algo-small" });
      const groups = pairs.map((pair, k) => {
        const x = left + 8 + k * step;
        pairText(svg, x + (n * (bar + 2)) / 2, base + 14, pair);
        const bars = data.profiles.map((_, s) => {
          const rect = el(svg, "rect", {
            x: x + s * (bar + 2), y: base - height, width: bar, height, class: "algo-bar",
          });
          rect.style.fill = sampleColor(s);
          return rect;
        });
        // The gap between the compared bars' tops, beside them
        const gap = el(svg, "rect", { x: x + n * (bar + 2) + 2, y: base, width: 5, height: 0, class: "algo-diff-fill" });
        return { bars, gap };
      });
      const grow = (fraction) => groups.forEach(({ bars }, k) => bars.forEach((rect, s) => {
        rect.style.transform = `scaleY(${fraction * data.profiles[s][k] / scale.top})`;
      }));

      // One row per comparison: the squared differences stacked into the
      // distance, which stays
      text(svg, distanceRows.left, 250, "Squared differences, stacked and square-rooted", {
        "text-anchor": "start", class: "algo-small",
      });
      const rows = distanceRows.draw(svg, 268);

      if (run.animate) {
        grow(0);
        await run.pause(50);
      }
      grow(1);
      await run.pause(900);
      for (const [r, [i, j]] of comparisons.entries()) {
        const row = rows[r];
        const color = pairColor(r);
        // The two samples compared stand out
        versus(heading, [i, j]);
        key.style.fill = color;
        keys.forEach((g, s) => g.classList.toggle("algo-dim", s !== i && s !== j));
        rows.forEach(({ label }, q) => label.classList.toggle("algo-dim", q !== r));
        groups.forEach(({ bars, gap }, k) => {
          bars.forEach((rect, s) => rect.classList.toggle("algo-dim", s !== i && s !== j));
          const a = data.profiles[i][k] / scale.top;
          const b = data.profiles[j][k] / scale.top;
          gap.setAttribute("y", base - Math.max(a, b) * height);
          gap.setAttribute("height", Math.abs(a - b) * height);
          gap.style.fill = color;
          gap.classList.add("algo-hidden");
        });
        await run.pause(900);
        for (const { gap } of groups) gap.classList.remove("algo-hidden");
        await run.pause(700);
        // Each difference leaves the chart, turning and stretching as it flies
        // into the row, where it becomes its share of the distance
        const stack = el(svg, "g");
        const segments = distanceRows.segments([i, j]).map(({ x, w }, k) => {
          const { gap } = groups[k];
          const x0 = Number(gap.getAttribute("x"));
          const y0 = Number(gap.getAttribute("y"));
          const w0 = Number(gap.getAttribute("width"));
          const h0 = Math.max(0.5, Number(gap.getAttribute("height")));
          const h = distanceRows.bar;
          const segment = el(stack, "rect", {
            x: x0, y: y0, width: w0, height: h0, class: "algo-diff-fill algo-slide",
          });
          segment.style.fill = color;
          // Turned a quarter, so the gap's height runs along the row
          segment.dataset.to = `translate(${x - x0}px, ${row.y - y0 + h}px) `
            + `rotate(-90deg) scale(${h / w0}, ${w / h0})`;
          return segment;
        });
        if (run.animate) await run.pause(50);
        for (const [k, segment] of segments.entries()) {
          groups[k].gap.classList.add("algo-hidden");
          segment.style.transform = segment.dataset.to;
          await run.pause(150);
        }
        await run.pause(1000);
        const d = data.distance[i][j];
        row.value.textContent = `d = √${(d ** 2).toFixed(4)} = ${d.toFixed(3)}`;
        await run.pause(2000);
        row.value.textContent = `d = ${d.toFixed(3)}`;
      }
      heading.replaceChildren();
      for (const node of svg.querySelectorAll(".algo-dim")) node.classList.remove("algo-dim");
    },

    // 5. The distances, as left by the comparison, become similarities
    async similarity(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 200");
      const { left, width } = distanceRows;
      const heading = text(svg, left, 14, "", { "text-anchor": "start" });
      const rows = distanceRows.draw(svg, 50).map((row, r) => {
        const stack = distanceRows.stack(svg, row, r);
        const rect = el(svg, "rect", {
          x: left, y: row.y, width, height: distanceRows.bar, class: "algo-bar algo-bar-x",
        });
        rect.style.fill = pairColor(r);
        rect.style.transform = "scaleX(0)";
        return { ...row, stack, rect };
      });
      const end = rows[rows.length - 1].y + distanceRows.bar + 6;
      const marker = el(svg, "line", { x1: left + width, x2: left + width, y1: 38, y2: end, class: "algo-marker" });
      const markerLabel = text(svg, left + width, 30, "max d", { class: "algo-small" });
      const cell = 44;
      const title = text(svg, 500 + (n * cell) / 2, 30, "Similarity");
      const setCell = drawMatrix(svg, 500, 60, cell);
      for (const node of [marker, markerLabel, title]) node.classList.add("algo-hidden");
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) setCell(i, j, 0, "");

      heading.textContent = "Distances";
      rows.forEach(({ value }, r) => {
        const [i, j] = comparisons[r];
        value.textContent = `d = ${data.distance[i][j].toFixed(3)}`;
      });
      await run.pause(1500);
      heading.textContent = "Divided by the largest distance";
      for (const node of [marker, markerLabel]) node.classList.remove("algo-hidden");
      rows.forEach(({ value }, r) => {
        const [i, j] = comparisons[r];
        value.textContent = `d / max d = ${(data.distance[i][j] / maxD).toFixed(2)}`;
      });
      await run.pause(1800);
      heading.textContent = "Subtracted from 1: similarities";
      markerLabel.textContent = "1";
      rows.forEach(({ rect, value, stack }, r) => {
        const [i, j] = comparisons[r];
        const s = data.similarity[i][j];
        stack.classList.add("algo-hidden");
        rect.style.transform = `scaleX(${s})`;
        value.textContent = `1 − d / max d = ${s.toFixed(2)}`;
      });
      await run.pause(1200);
      title.classList.remove("algo-hidden");
      for (let i = 0; i < n; i++) setCell(i, i, 1, "1");
      for (const [i, j] of comparisons) {
        const s = data.similarity[i][j];
        setCell(i, j, s, s.toFixed(2));
        setCell(j, i, s, s.toFixed(2));
        await run.pause(300);
      }
    },

    // only_transitions: links within a label are dropped
    async transitions(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 260");
      const sample = data.samples[0];
      const scale = fit(data.samples, 360, 240);
      const { links } = drawSample(svg, sample, [0, 10, 380, 240], scale, true);
      const max = Math.max(...sample.counts.flat());
      const setCell = drawTable(svg, 450, 80, 50, `${data.names[0]}: connection counts`);
      const show = (dropSame) => labels.forEach((_, i) => labels.forEach((_, j) => {
        const v = dropSame && i === j ? 0 : sample.counts[i][j];
        setCell(i, j, v / max, String(v));
      }));
      show(false);
      for (const link of links) link.classList.add("algo-counted");
      await run.pause(1200);
      sample.pairs.forEach(([i, j], k) => {
        if (sample.label[i] === sample.label[j]) links[k].classList.add("algo-hidden");
      });
      show(true);
    },

    // pools: two samples' tables are summed
    async pools(svg, run) {
      svg.setAttribute("viewBox", "0 0 640 220");
      const [a, b] = data.samples;
      const total = a.counts.map((row, i) => row.map((v, j) => v + b.counts[i][j]));
      const max = Math.max(...total.flat());
      const cell = 40;
      for (const [counts, name, x] of [
        [a.counts, data.names[0], 45],
        [b.counts, data.names[1], 255],
      ]) {
        const setCell = drawTable(svg, x, 80, cell, name);
        counts.forEach((row, i) => row.forEach((v, j) => setCell(i, j, v / max, String(v))));
      }
      text(svg, 205, 80 + 1.5 * cell, "+", { class: "algo-big" });
      const equals = text(svg, 415, 80 + 1.5 * cell, "=", { class: "algo-big" });
      const setPool = drawTable(svg, 465, 80, cell, "Pool");
      equals.classList.add("algo-hidden");
      labels.forEach((_, i) => labels.forEach((_, j) => setPool(i, j, 0, "")));
      await run.pause(1000);
      equals.classList.remove("algo-hidden");
      for (let i = 0; i < labels.length; i++) {
        for (let j = 0; j < labels.length; j++) {
          setPool(i, j, total[i][j] / max, String(total[i][j]));
          await run.pause(150);
        }
      }
    },
  };

  const descriptions = {
    spots: "Three tissue sections of labelled spots.",
    links: "Links drawn between each spot and its neighbors.",
    counts: "Each sample's links counted in turn into its label by label table.",
    profile: "Each sample's table as bars of counts, then of proportions, on shared axes.",
    compare: "All samples' profiles side by side, with each pair's squared differences stacked in turn into one bar.",
    similarity: "The distances, divided by the largest and subtracted from 1, filling a similarity matrix.",
    transitions: "Links within a label dropped, zeroing the table's diagonal.",
    pools: "Two samples' tables summed into one pool.",
  };

  for (const container of document.querySelectorAll(".algorithm-figure")) {
    const name = container.dataset.figure;
    const draw = figures[name];
    if (!draw) continue;
    const svg = el(container, "svg", { role: "img", "aria-label": descriptions[name] });
    let run = null;
    const play = (animate) => {
      if (run) run.cancelled = true;
      run = new Run(animate);
      svg.replaceChildren();
      svg.classList.remove("algo-pending");
      svg.classList.toggle("algo-still", !animate);
      draw(svg, run, container.dataset).catch((e) => {
        if (e !== CANCELLED) throw e;
      });
    };
    if (still || !("IntersectionObserver" in window)) {
      play(false);
      continue;
    }
    const button = document.createElement("button");
    button.type = "button";
    button.className = "algo-replay";
    button.textContent = "Replay";
    button.addEventListener("click", () => play(true));
    container.append(button);
    // Drawn but hidden until first scrolled into view (or its tab opened),
    // so it keeps its size
    play(false);
    svg.classList.add("algo-pending");
    const observer = new IntersectionObserver((entries) => {
      if (entries.some((e) => e.isIntersecting)) {
        observer.disconnect();
        play(true);
      }
    }, { threshold: 0.5 });
    observer.observe(container);
  }
})();
