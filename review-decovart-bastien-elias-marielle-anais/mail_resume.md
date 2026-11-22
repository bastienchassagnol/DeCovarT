Objet : DeCovarT, réponses à tes remarques

Bonjour Elias (Marielle, Anaïs en copie),

Merci pour ta relecture. Je te joins une note (PDF et HTML) qui reprend tes remarques une par une, avec les liens vers les vignettes et le livre de reproductibilité. En bleu tes remarques, en vert ce que je reprends, en orange mes nuances. Les points principaux :

1. **Gènes nuls dans un type.** D'accord pour la prédiction, mais je garde la porte de présence. Avec un facteur de plateforme par gène (ton biais multiplicatif), des marqueurs parfaits rendent la composition non identifiable (contre-exemple de MEAD). Une variance nulle rend la covariance du type singulière, mais pas celle du mélange à l'intérieur du simplexe.

2. **Covariance et pondération.** Tu as raison. Le score se coupe en un quasi-score GLS et un score de covariance, et la part d'information de ce dernier décroît en 1/n : à grand n, DeCovarT équivaut à une version « covariance pleine » de MuSiC. Deux réserves : avec une couche entre réplicats (ta formule, qui est notre couche 4), cette part ne s'annule plus et n_i devient un nombre effectif borné ; sur Suppinger, on estime n̂ < 1 cellule à 96 h.

3. **Régularisation et critère.** Le seuil G+1 est utile, d²log G beaucoup moins, puisque le plug-in (shrinkage de Schäfer–Strimmer) ne met aucune corrélation à zéro. Le niveau choisi par la vraisemblance des bulks (0,1 à 0,5) est bien plus faible que sur les cellules (0,66 à 0,94), sans jamais atteindre la diagonale. Je propose des arguments pour défendre la vraisemblance contre AIC/BIC au sein d'un bras ; c'est une question ouverte.

4. **Design.** En relisant Suppinger, les bulks viennent d'une troisième lignée (129sv/ev), différente des deux lignées single-cell, et la référence n'a pas de réplicats biologiques. Le jeu de Villaronga-Luque et al. (2025), avec des structures individuelles marquées par MULTI-seq, permettrait d'estimer la couche entre gastruloïdes.

5. **Simulations.** Le scénario single-cell actuel n'a ni biais par gène, ni effet donneur, ni comptages. Je propose de les ajouter et de faire de f_cov (la part d'information portée par la covariance) un facteur explicite, plutôt que le seul recouvrement. On comparerait quatre estimateurs emboîtés : moyenne seule, GLS diagonal, GLS plein et DeCovarT.

Trois points me tiennent particulièrement à cœur : le critère de sélection, la preuve asymptotique pour le modèle single-cell, et ce plan de simulation. Ton avis m'aiderait beaucoup.

Bien à toi,
Bastien
