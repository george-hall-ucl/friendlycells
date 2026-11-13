// Shows the R or Python wording of text marked {.lang-r} or {.lang-py},
// following the language tabs. Quarto restores the last chosen tab from
// localStorage without firing events, so read the same setting here.
(() => {
  const show = (tab) => {
    document.documentElement.dataset.lang = tab === "Python" ? "py" : "r";
  };
  let saved = "R";
  try {
    const tabs = localStorage.getItem("quarto-persistent-tabsets-data");
    saved = JSON.parse(tabs || "{}").language || "R";
  } catch (e) {
    // Storage blocked: fall back to R, the first tab
  }
  show(saved);
  document.addEventListener("click", (e) => {
    const tab = e.target.closest('div[data-group="language"] a[id^="tabset-"]');
    if (tab) show(tab.textContent.trim());
  });
})();
