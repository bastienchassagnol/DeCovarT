
Bonjour Elias,




Merci pour ta relecture détaillée et encore désolé pour le délai. Tes remarques fort pertinentes, dignes d'un "relecteur 2" n'appelaient pas de réponses simples.

Je joins par ailleurs une note PDF qui reprend chaque point

(en bleu tes remarques, en vert mes accords, en orange mes nuances) ; les sections à lire sont indiquées entre crochets,

voire aussi sa version versionnée: https://github.com/bastienchassagnol/DeCovarT/tree/main/review-decovart-bastien-elias-marielle-anais.





## 1. Ce que tu as relevé





Points positifs : le cadrage en erreurs sur les variables, le choix d'une covariance robuste plutôt que d'un réseau,

et l'idée que l'apport principal d'un modèle multivarié est une incertitude bien calibrée.





Points problématiques :





- **Gènes marqueurs.** Faut-il vraiment exclure les gènes exprimés dans un seul type cellulaire (ma « presence gate ») ? [§2.1 « Faut-il garder les gènes marqueurs parfaits ? »]

- **Régularisation de la covariance.** Comment choisir le niveau de régularisation λ

de la covariance de chaque type sans sur-apprentissage, puis comparer des modèles différents ?

[§2.3 « Choisir λ, puis comparer des modèles »]

- **Apport réel de la covariance.** À grand nombre de cellules, la covariance n'aiderait

que par la pondération des gènes, comme une version « covariance pleine » de MuSiC. [§3.3

« La log-vraisemblance : une distance de Mahalanobis plus un déterminant »]

- **Variance cachée entre réplicats.** Une variabilité entre gastruloïdes, absente du modèle,

deviendrait certainement le terme de variabilité dominant (se rapprochant ainsi de la première version de Decovart).

[§3.5 « Applications sur des données réelles »]

- **Simulations trop optimistes** sans biais par gène, effet donneur ni couche de comptage.

[§3.4 « Applications sur des données synthétiques »]





## 2. Mes réponses provisoires, à discuter





1. **Marqueurs : pas de réponse universelle.** Avec DeCovarT (une seule échelle par échantillon),

les marqueurs sont identifiables et très informatifs.

Avec un facteur libre par gène, comme dans MEAD, ils rendent le modèle non identifiable

à partir des moyennes seules. Les garder pose trois difficultés pratiques :

la vraisemblance n'est plus bornée sur les bords du simplexe (type absent ou dominant),

la covariance par type devient singulière, et le biais de plateforme se concentre sur eux.

Je garde pour l'instant la presence gate, et je propose de comparer les performances avec et sans marqueurs,

avec un plancher de variance poissonien. [§2.1, arbre de décision]





2. **Régularisation : aucune des trois stratégies n'est parfaite.** J'ai comparé une validation croisée

sur les cellules, des pseudo-bulks et la vraisemblance des vrais bulks.

Je retiens la dernière (λ entre 0,1 et 0,5), en sachant qu'elle est optimiste :

λ est choisi sur les bulks qu'il sert à déconvoluer.





Les deux pistes plus propres demandent soit beaucoup plus de réplicats biologiques à composition voisine

(validation croisée sur les bulks type LOOCV), soit de traiter λ comme un hyperparamètre estimé dans le modèle.

Pour la seconde, on pourrait envisager une approche bayésienne complète (un a priori sur λ) ou

la vraisemblance marginale en λ (Bayes empirique).









2bis: Pour comparer la performance des modèles entre eux, en plus de la vraisemblance,

l'UACV me semble plus réaliste que l'AIC ou le BIC, avec seulement trois réplicats. [§2.3]





3. **Apport de la covariance : tu as raison asymptotiquement.** La part d'information portée par la covariance décroît en 1/n_i,

avec n_i le nombre de **cellules** du bulk. À grand n_i, DeCovarT équivaut donc à

un GLS à covariance pleine, MuSiC n'en étant que la version diagonale (WLS).





Deux remarques : sur Suppinger, n_i estimé est inférieur à une cellule, loin de ce régime ; et ce second terme

ne tend plus vers zéro dès lors qu'on ajoute une variance entre gastruloïdes. [§3.3, §3.2]





4. **Variance entre gastruloïdes : à modéliser, mais pas avec Suppinger.** Ce terme est quadratique en p_j ;

avec 10⁴ à 10⁶ cellules il domine le terme linéaire dès qu'il n'est pas négligeable.





Malheureusement, la référence single-cell de Suppinger

n'a qu'un réplicat biologique.

Le jeu de Villaronga-Luque et al. (2025), avec des gastruloïdes individuels marqués par MULTI-seq,

permettrait d'estimer cette variance, et de tester le modèle à deux niveaux

en comparant les analyses compositionnelles estimées avec celles effectuées sur les données

de scRNASeq directement (scCODA, sccomp). [§3.5, encadré « Sources de variabilité biologique »]





5. **Simulations : garder le plan actuel, mais rapporter f_cov.** Je garde les scénarios actuels,

où le critère d'overlap est la principale variable d'ajustement(un bon indicateur global de difficulté,

mais qui ne distingue pas les méthodes). Toutefois, je rapporte en plus pour chaque scénario le ratio de la part d'information portée par la covariance, f_cov, liée

au choix d'une loi de distribution gaussienne multivariée, sur f_exp (distance de Mahalanobis).





Un plan qui fait varier f_cov à information moyenne fixée n'a, à ma connaissance, jamais été proposé :

ce serait une contribution méthodologique en soi, à garder pour plus tard.





Enfin, scDesign3 servirait ensuite à tester la robustesse en utilisant d'autres lois de comptage. [§3.4, annexe 4.3]





J'écarterais en revanche pour le moment l'atlas d'embryons (Pijuan-Sala et al.) pour tester la robustesse du non-alignement

entre donneurs de bulk et fournisseurs des profils de single-cell: la biologie et le calendrier de développement diffèrent

de ceux des gastruloïdes, et ce désaccord entre référence et bulk brouillerait la lecture des résultats.

En effet, l'ajout requis d'une covariable temporelle (requise si développement asynchrone) pour aligner les deux temps de développement

pourrait tout aussi bien capturer l'effet batch.





## Conclusion





Je pense qu'on a déjà bien souligné les limites des méthodes actuelles, et qu'on propose une contribution certes imparfaite, mais qui est le premier

modèle à proposer un vrai test flexible pour comparer l'impact de covariables (condition biologique, temps, ...) sur al composition cellulaire

Complexifier encore le modèle ou ajouter beaucoup de simulations risquerait selon moi de desservir l'article, en fonction du

niveau d'exigence dépend du journal visé.





(disclaimer: demandé par IA pour les journaux suivants)

- **Bioinformatics** : format court (article d'environ sept pages) et lectorat méthodologique ;

publication la plus rapide et la moins exigeante en validation biologique,

mais visibilité moindre et peu de place pour la théorie, qui partirait en supplément. Bien aligné sur le package et le benchmark, moins sur les apports théoriques.

- **PLOS Computational Biology** (la cible de CardamomOT) : pas de limite de longueur, mais frais de publication élevés.

Place pour la théorie et le logiciel dans le corps du texte ; impact intermédiaire.

Très bien aligné sur DeCovarT tel qu'il est.

- **Nature Methods** : visibilité et impact les plus forts, mais format très contraint

(environ 3 000 mots, six figures ou tableaux, cinquante références) + exigence d'une nouveauté large et

d'une validation sur données réelles avec vérité de terrain, que Suppinger ne fournit pas ; la théorie irait en supplément.

Il est par ailleurs fort possible que les limites soulevées précédemment soient rédhibitoires si identifiées par les relecteurs.





Vos avis sur ces cinq points, et sur le choix du journal, m'aideraient beaucoup 😊
