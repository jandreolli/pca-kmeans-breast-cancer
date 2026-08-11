```r
# TP Noté Fondements Statistiques - ANDREOLLI Justine M1 SDSC
#
# Analyse du jeu de données METABRIC :
# - Préparation des données
# - ACP et k-means
# - Interprétation des clusters
# - Prédiction de l'indice de Nottingham


# ============================================================
# 1. PRÉPARATION DU JEU DE DONNÉES
# ============================================================

library(tidyverse)

metabric.data <- read.csv2(
  './data/METABRIC_RNA_Mutation.csv',
  header = TRUE,
  sep = ",",
  row.names = 'patient_id'
)

# Supprimer les données de mutation que nous n'allons pas utiliser
genes.data <- metabric.data[, 1:519]

# Conversion des colonnes au bon type
genes.data <- type_convert(genes.data)

# Première diminution de dimension par sélection des 100 gènes
# les plus corrélés au score de gravité de Nottingham

corre <- c()

for (i in 1:489) {
  corre <- c(
    corre,
    cor(genes.data[, 21], genes.data[, i + 30])^2
  )
}

selected <- order(corre, decreasing = TRUE)[1:100]

genes.data <- genes.data[, c(1:30, 30 + selected)]


# ============================================================
# 2. ACP ET K-MEANS
# ============================================================

# ------------------------------------------------------------
# Question 1
# ------------------------------------------------------------
# Réaliser une ACP normée sur les données des 100 gènes en
# gardant les variables médicales comme variables
# supplémentaires d'interprétation.
#
# Les variables 1, 19, 20, 21, 23 et 28 sont quantitatives.
# Les autres variables sont qualitatives.

library(FactoMineR)
library(factoextra)

acp_rslt <- PCA(
  genes.data[, -(1:30)],
  scale.unit = TRUE,
  quanti.sup = c(1, 19, 20, 21, 23, 28),
  quali.sup = setdiff(1:30, c(1, 19, 20, 21, 23, 28))
)

fviz_screeplot(
  acp_rslt,
  addlabels = TRUE,
  ylim = c(0, 50)
)

# Variance expliquée par chaque composante
variance_expliquee <- acp_rslt$eig[, 2]

# Variance cumulée
variance_cumulee <- acp_rslt$eig[, 3]

# Nombre de composantes permettant d'atteindre au moins 80 %
nb_composantes <- which(variance_cumulee >= 80)[1]

cat(
  "\nNombre de composantes gardées :",
  nb_composantes,
  "\n"
)

# D'après l'analyse, 34 composantes sont nécessaires pour
# atteindre environ 80.1 % de la variance totale.
#
# Cela s'explique par le fait que chaque composante possède
# une variance individuelle relativement faible.
# La première composante explique environ 19.3 % de la
# variance, puis les suivantes contribuent progressivement
# moins à la variance totale.


# ------------------------------------------------------------
# Question 2
# ------------------------------------------------------------
# Appliquer l'algorithme des k-means sur les vecteurs projetés,
# pour un nombre de composantes allant de 1 à 30.
#
# On étudie l'évolution de la somme des carrés intra-classes
# afin de choisir le nombre de clusters avec le critère du coude.

library(cluster)
library(factoextra)

coord_projetees <- acp_rslt$ind$coord

# Nombre de dimensions disponibles
nb_dimensions <- ncol(coord_projetees)

# Somme des carrés intra-clusters
wss <- c()

for (k in 1:30) {
  
  kmeans_rslt <- kmeans(
    coord_projetees[, 1:min(30, nb_dimensions)],
    centers = k,
    nstart = 10,
    iter.max = 100
  )
  
  wss <- c(wss, kmeans_rslt$tot.withinss)
}

# Visualisation du critère du coude
plot(
  1:30,
  wss,
  type = "b",
  pch = 19,
  frame = FALSE,
  xlab = "Nb de clusters",
  ylab = "Within-Cluster Sum of Squares (WSS)",
  main = "Critère du coude"
)

# Choix de 3 clusters à partir du critère du coude
nombre_optimal_clusters <- 3

# Application finale du k-means avec 3 clusters
final_kmeans <- kmeans(
  coord_projetees[, 1:min(30, nb_dimensions)],
  centers = nombre_optimal_clusters,
  nstart = 10,
  iter.max = 100
)

# Ajouter le cluster de chaque patiente au jeu de données
genes.data$cluster <- as.factor(final_kmeans$cluster)

# Avec l'algorithme k-means et 3 clusters,
# les patientes sont classées en 3 groupes.


# ------------------------------------------------------------
# Question 3
# ------------------------------------------------------------
# Tracer le nuage de points dans le plan principal de l'ACP
# avec des couleurs correspondant aux classes de k-means.

fviz_pca_ind(
  acp_rslt,
  geom = "point",
  col.ind = genes.data$cluster,
  palette = c("#2E9FDF", "#996DCF", "#850606"),
  addEllipses = TRUE,
  ellipse.type = "convex",
  legend.title = "Cluster",
  repel = TRUE
) +
  ggtitle(
    "Nuage de points dans le plan principal de l'ACP\navec les clusters de k-means"
  )

# On utilise les deux premières composantes pour visualiser
# les dimensions.
#
# On remarque la présence de 3 clusters.
# Les clusters 1 et 2 se confondent légèrement dans certaines
# zones. Cela pourrait signifier que d'autres composantes
# permettraient de mieux les séparer.
#
# La dimension 1 est plus discriminante que la dimension 2 :
# environ 19.3 % de variance expliquée pour la première
# composante contre environ 6.5 % pour la deuxième.


# ============================================================
# 3. INTERPRÉTATION DES RÉSULTATS
# ============================================================

# ------------------------------------------------------------
# Question 4
# ------------------------------------------------------------
# À l'aide de tests du chi-deux, rechercher les variables
# médicales qualitatives auxquelles les classes trouvées sont
# les plus liées.
#
# On compare ensuite les trois variables les plus significatives
# avec la classification obtenue à l'aide de tables de
# contingence.

variables_qualitatives <- setdiff(
  1:30,
  c(1, 19, 20, 21, 23, 28)
)

rslt_chi2 <- data.frame(
  Variable = character(),
  p_value = numeric()
)

# Test du Chi-deux pour chaque variable qualitative
for (var in variables_qualitatives) {
  
  table_contingence <- table(
    genes.data[[var]],
    genes.data$cluster
  )
  
  test_chi2 <- chisq.test(table_contingence)
  
  rslt_chi2 <- rbind(
    rslt_chi2,
    data.frame(
      Variable = colnames(genes.data)[var],
      p_value = test_chi2$p.value
    )
  )
}

# Trier les variables par p-value croissante
rslt_chi2 <- rslt_chi2[
  order(rslt_chi2$p_value),
]

print(rslt_chi2)

# Comparaison des trois premières variables qualitatives
for (var in rslt_chi2$Variable[1:3]) {
  
  cat(
    "\nTable de contingence de la variable :",
    var,
    "\n"
  )
  
  print(
    table(
      genes.data[[var]],
      genes.data$cluster
    )
  )
}

# Les variables les plus significatives, avec les p-values
# les plus faibles, sont :
# - integrative_cluster
# - pam50_._claudin.low_subtype
# - X3.gene_classifier_subtype
#
# Pour la table de contingence de integrative_cluster,
# la modalité 10 est un bon indicateur du cluster 3.
#
# Pour pam50_._claudin.low_subtype, la modalité Basal est
# fortement liée au cluster 3 et LumA au cluster 2.
#
# Pour X3.gene_classifier_subtype, la modalité ER-/HER2-
# est principalement associée au cluster 3 et
# ER+/HER2- Low Prolif au cluster 2.


# ============================================================
# 4. PRÉDICTION DE L'INDICE DE NOTTINGHAM
# ============================================================

# ------------------------------------------------------------
# Question 5
# ------------------------------------------------------------
# Mettre en place une prédiction de l'indice de Nottingham
# à l'aide d'un modèle linéaire.

# Variable cible
y <- genes.data$nottingham_prognostic_index

# Variables explicatives :
# - le cluster obtenu avec k-means
# - les 100 gènes sélectionnés
X <- genes.data[
  ,
  c("cluster", colnames(genes.data)[31:130])
]

X$cluster <- as.factor(X$cluster)

# Séparation des données en ensemble d'entraînement et de test
set.seed(48)

train_index <- sample(
  1:nrow(X),
  0.8 * nrow(X)
)

X_train <- X[train_index, ]
y_train <- y[train_index]

X_test <- X[-train_index, ]
y_test <- y[-train_index]

# Application du modèle linéaire
model <- lm(
  y_train ~ .,
  data = X_train
)

summary(model)

# Prédictions sur l'ensemble de test
y_pred <- predict(
  model,
  newdata = X_test
)

# Calcul du RMSE
rmse <- sqrt(
  mean((y_test - y_pred)^2)
)

cat(
  "RMSE sur l'ensemble de test :",
  rmse,
  "\n"
)

# L'indice de Nottingham étant une variable continue,
# la régression linéaire est utilisée.
#
# Les clusters de k-means sont inclus comme variables
# explicatives car ils regroupent les patientes selon
# des caractéristiques similaires, ce qui peut contribuer
# à expliquer les variations de l'indice de Nottingham.
#
# Environ 25 % de la variance est expliquée avec une erreur
# quadratique moyenne d'environ 0.998.
#
# Les clusters de k-means et les gènes sélectionnés permettent
# de repérer plusieurs variables importantes, notamment :
# - gsk3b
# - ncoa3
# - mmp12
# - tp53bp1
#
# Cependant, des modèles non linéaires, comme Random Forest,
# pourraient être utilisés pour améliorer les performances,
# puisque le modèle actuel n'explique qu'environ un quart
# de la variance totale de l'indice de Nottingham.
```
