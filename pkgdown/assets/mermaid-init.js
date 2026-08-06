/**
 * pkgdown Mermaid renderer (shared after_body include).
 * Handles pkgdown/README <pre class="mermaid"> and Quarto
 * <pre class="mermaid mermaid-js"> / language-mermaid fences.
 */
(function () {
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
    return (code ? code.textContent : pre.textContent) || "";
  }

  async function run() {
    if (!window.mermaid) return;
    var pres = collect();
    if (!pres.length) return;
    // Core mermaid themes: default | dark | forest | neutral | base
    // (Quarto's "sandstone" is not a mermaid.js theme.)
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
