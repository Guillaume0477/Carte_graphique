# Carte_graphique

# TP 3 relatif au cartes graphiques et compute shaders
# Auteurs : DURET Guillaume

Le ficher tutos/M2/tuto_raytrace_compute.cpp represente le code final
(rendu pour le 31 janvier) du TP3 du module de programmation sur cartes graphiques.


update 1/02
ajout du bvh cousu et parcourt avec seulement skip (il est possible d'avoir next en faisant -1 à l'indice du noeud courant)
mais les feuilles utilisent toujours left et right...
Les performances sont dans le meme ordre de grandeur que le parcourt du bvh avec une faible pile ??
