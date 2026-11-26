# Retour figures 1 et 4 (Elias, Marielle, Anaïs)

- **HTML:** open `retour_figure_1_and_4.html` in a browser (needs the
  `retour_figure_1_and_4_files/` folder beside it). Figures 1 et 4 du
  manuscrit ne sont pas embarquées : liens vers
  `../article/figs/fig01_pipeline.pdf` et
  `../article/figs/fig04_suppinger.pdf`.
- **PDF:** `retour_figure_1_and_4.pdf` est la version Ghostscript
  (`/ebook`, 150 dpi). L’export LuaLaTeX brut est gardé sous
  `retour_figure_1_and_4_render.pdf` (gitignored).

```bash
cd review-decovart-bastien-elias-marielle-anais
./render_retour_figure_1_and_4.sh
# plus léger : PDFSETTINGS=/screen ./render_retour_figure_1_and_4.sh
```

Compression manuelle (ne pas écraser l’entrée) :

```bash
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/ebook \
  -dNOPAUSE -dQUIET -dBATCH \
  -sOutputFile=retour_figure_1_and_4.pdf retour_figure_1_and_4_render.pdf
```
