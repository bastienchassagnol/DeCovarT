/**
 * pkgdown Mermaid renderer for the home page (README.md fences).
 * Vignettes are Quarto-rendered to SVG at build time; index.html still
 * needs a client-side pass for <pre class="mermaid"> blocks.
 */
(function () {
  function stripQuartoInit(text) {
    // Quarto themes such as "sandstone" are not mermaid.js themes.
    return String(text || "")
      .replace(/^\s*%%\{init:[\s\S]*?\}%%\s*/m, "")
      .trim();
  }

  function collect() {
    var out = [];
    document.querySelectorAll("pre.mermaid, pre.mermaid-js").forEach(function (pre) {
      if (pre.querySelector("svg")) return;
      out.push(pre);
    });
    document.querySelectorAll(
      "pre code.language-mermaid, pre code.mermaid, pre code[class*='language-mermaid']"
    ).forEach(function (code) {
      var pre = code.closest("pre");
      if (!pre || pre.querySelector("svg")) return;
      if (out.indexOf(pre) === -1) out.push(pre);
    });
    return out;
  }

  function textOf(pre) {
    var code = pre.querySelector("code");
    return stripQuartoInit(code ? code.textContent : pre.textContent);
  }

  async function run() {
    if (!window.mermaid) return;
    var pres = collect();
    if (!pres.length) return;
    window.mermaid.initialize({
      startOnLoad: false,
      theme: "neutral",
      securityLevel: "loose",
      flowchart: { htmlLabels: true, curve: "basis" }
    });
    var nodes = pres.map(function (pre, i) {
      var div = document.createElement("div");
      div.className = "mermaid";
      div.id = "pkgdown-mermaid-" + i;
      div.textContent = textOf(pre);
      pre.replaceWith(div);
      return div;
    });
    try {
      await window.mermaid.run({ nodes: nodes });
    } catch (e) {
      console.warn("Mermaid render failed:", e);
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", run);
  } else {
    run();
  }
})();
