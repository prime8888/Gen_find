# Projet Gen Find

## Membres
- Daniil Kudriashov
- Reem Baza
- Kahina Nessah

## Description du projet
Le projet Gen Find vise à extraire les régions fonctionnelles des génomes de différents organismes récupérés depuis la base de données genome de la National Library of Medicine. Les séquences récupérées sont stockées localement dans des fichiers texte par organisme, région fonctionnelle et par ID de record. Le programme utilise le multithreading, ce qui permet d'améliorer la performance et d'optimiser les ressources du système. Cependant, comme les données de génomes sont lourdes, leur téléchargement prend beaucoup de temps. De plus, soyez conscient que les fichiers générés peuvent prendre beaucoup d'espace disque (>50 Go par famille d'organismes).

## Fonctionnalités
Les fonctionnalités implémentées par rapport aux requis dans le sujet du projet incluent :
- Choix des régions fonctionnelles : CDS, centromère, intron, mobile_element, ncRNA, rRNA, télomère, tRNA, 3'UTR, 5'UTR
- Menu permettant le choix d'un ensemble d'organismes
- Interface graphique
- Génération et mise à jour de l'arborescence locale
- Tests sur les bornes
- Opérateurs : join, complement, complement(join)
- Fenêtre log dans l'interface avec des messages d'information
- Un seul bouton d'exécution pour gérer l'arrêt et le démarrage du parsing
- Interface pouvant être déplacée avec une taille modifiable
- Le programme ne se bloque pas sur les données corrompues
- Les organismes déjà parsés ne sont pas reparsés lors des exécutions suivantes

## Lancer le projet à partir des fichiers source

### Requirements
- Python 3.11
- Les librairies nécessaires peuvent être trouvées dans le fichier `environment.yml`

### Exécution
Après avoir installé les dépendances, exécutez le fichier `main.py` se trouvant dans le répertoire `src/`.

## Structure des fichiers
- `src/` - fichiers source du code
- `Results/` - arborescence avec les résultats du parsing (générée automatiquement lors du premier lancement)
- `local_overview.txt` - fichier récupéré depuis NCBI avec les organismes disponibles pour construire l'arborescence en local (généré automatiquement lors du premier lancement)
- `README.md` - ce fichier

## Exécutable
L'exécutable peut être téléchargé au lien suivant (version Windows) :
[Gen Find Windows Version](https://drive.google.com/file/d/1lnUwcTYWRMtD3B4bYQMYq2UpsrZm6fjM/view?usp=sharing)
NB : L'antivirus peut interdire le fichier, il faut donc l'ajouter aux exclusions.

## Utilisation du programme
Après le premier lancement, l'arborescence sera générée. Vous pouvez ensuite choisir les organismes ou les familles d'organismes à parser. Il est également possible de choisir plusieurs organismes à la fois avec `Ctrl` ou `Shift`, comme dans votre explorateur de fichiers. Vous devez ensuite choisir au moins une région fonctionnelle, puis appuyer sur le bouton "Start Parsing". Vous verrez le progrès dans la section information et les logs dans la fenêtre des logs. Si vous souhaitez arrêter le processus, vous pouvez appuyer sur le bouton "Stop Parsing". Le programme attendra la fin des threads déjà en cours de travail, mais ne commencera pas de nouveaux records. À la fin de l'exécution, vous pouvez retrouver les résultats dans le répertoire `Results/` dans les dossiers respectifs des organismes.
