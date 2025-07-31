# Mymatlab1

## 🎯 Description

Ce dépôt regroupe un **pipeline MATLAB complet** pour reconstruire la matrice de densité réduite d’ état photoélectronique à partir de données de spectroscopie LDPI (Laser‑Driven PhotoElectron Imaging). Il inclut la simulation, la modulation IR, l’inversion régularisée (Tikhonov + TV), le bootstrap d’incertitude, et la visualisation des résultats (module/phase, écart‑type).

---

## Contenu du dépôt

- `script_principal.m` *(à renommer si besoin)* : pipeline complet – de la génération de densité à l’analyse finale.
- `RhoDE_2_RhoE.m` et `RhoE_2_RhoDE.m` : fonctions de conversion entre les bases énergie‑délai et énergie pure.
- `Sampling_Sami.m` : fonction de calcul de la taille d’échantillonnage (`NE`) selon la méthode Sami.
- `constants.mat` : fichier contenant toutes les constantes physiques (eV, fs, c, TW…).
- Autres scripts annexes :
  - `Phase_modulator.m` : calcul de la phase d’accumulation IR.
  - `eigen_decomposition.m`  
  - `script_dependantp.m` *(non écrasé lors des mises à jour automatisées)*

---

## Prérequis

- MATLAB (version ≥ R2019a) avec le toolbox Parallel Computing (optionnelle).
- Tous les fichiers `.m` et `constants.mat` dans le même dossier.
- Une installation Git opérationnelle pour le versioning.

---

## 🚀 Instructions d’utilisation

1. **Cloner le dépôt** :
   ```bash
   git clone https://github.com/SamiBoujra/Mymatlab1.git
   cd Mymatlab1
