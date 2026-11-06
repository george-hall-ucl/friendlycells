// Shows the code snippet for the chosen language and data. The language
// buttons click the page's R/Python tabs, so lang.js and the tabs stay in sync.
(() => {
  const selector = document.querySelector(".selector");
  if (!selector) return;
  const root = document.documentElement;
  const TAB = 'div[data-group="language"] a[id^="tabset-"]';

  const mark = () => {
    for (const button of selector.querySelectorAll("button")) {
      const choice = button.parentElement.dataset.choice;
      const on = button.dataset.value === root.dataset[choice];
      button.classList.toggle("active", on);
      button.setAttribute("aria-pressed", on);
    }
  };

  const setTech = (tech) => {
    root.dataset.tech = tech;
    try {
      localStorage.setItem("friendlycells-tech", tech);
    } catch (e) {
      // Storage blocked: the choice lasts until the page is left
    }
    mark();
  };

  const setLang = (lang) => {
    const label = lang === "py" ? "Python" : "R";
    const tabs = [...document.querySelectorAll(TAB)].filter(
      (a) => a.textContent.trim() === label,
    );
    // A tab below the selector, so Quarto keeps the selector in place
    const tab =
      tabs.find(
        (a) => selector.compareDocumentPosition(a) & Node.DOCUMENT_POSITION_FOLLOWING,
      ) || tabs[0];
    if (tab) tab.click();
    else root.dataset.lang = lang;
  };

  selector.addEventListener("click", (e) => {
    const button = e.target.closest("button");
    if (!button) return;
    if (button.parentElement.dataset.choice === "lang") setLang(button.dataset.value);
    else setTech(button.dataset.value);
  });
  new MutationObserver(mark).observe(root, { attributeFilter: ["data-lang"] });

  let saved = null;
  try {
    saved = localStorage.getItem("friendlycells-tech");
  } catch (e) {
    // Storage blocked: start with Visium
  }
  const techs = [...selector.querySelectorAll('[data-choice="tech"] button')].map(
    (b) => b.dataset.value,
  );
  setTech(techs.includes(saved) ? saved : techs[0]);
})();
