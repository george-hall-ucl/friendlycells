// Links function names inside code blocks to their section of the Reference.
// The code is highlighted at build time, so this is done in the page.
(() => {
  const names = ["spatial_similarity", "plot_neighborhoods"];
  const pattern = new RegExp(`\\b(${names.join("|")})\\b`);
  for (const code of document.querySelectorAll("pre code")) {
    const walker = document.createTreeWalker(code, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    for (const node of nodes) {
      // With a capturing group, odd parts are the names
      const parts = node.textContent.split(pattern);
      if (parts.length === 1 || node.parentElement.closest("a")) continue;
      const linked = parts.map((part, i) => {
        if (i % 2 === 0) return part;
        const a = document.createElement("a");
        a.href = `reference.html#${part}`;
        a.textContent = part;
        return a;
      });
      node.replaceWith(...linked);
    }
  }
})();
