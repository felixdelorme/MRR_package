# Projet : MK Packaged v2

## Description du Projet
Ce projet traite des données liées au Micro Rain Radar, mais aussi au croisement de données avec des données méteo (precipitation, temperatures, albedo). Il génère diverses figures et images, traite les données et fournit des visualisations au format Excel.

## Structure du Répertoire
  ```bash
    MRR_package/
    ├── Dat_processing/
    │   ├── Data/
    │   ├── albedo_figures/
    │   ├── precip_albedo_figures/
    │   ├── precip_temp_figures/
    │   ├── temperature_figures/
    │   ├── albedo.csv
    │   ├── Read_Dat.py
    │   ├── TestRead.py
    │
    ├── Melting_Layer/
    │   ├── MRR_ML_detection_auto.py
    │   ├── MRR_ML_detection.py
    │   ├── MRR_ploting.py
    │
    ├── Micro_Rain_Radar/
    │   ├── .ipynb_checkpoints/
    │   ├── Data/
    │   ├── lib/
    │   ├── MK_processed/
    │   │   ├── 202309/
    │   │   ├── 202310/
    │   │   ├── 202311/
    │   │   ├── PYRAMIDE/
    │   ├── default_parameters.txt
    │   ├── Run_MK.py
    │
    ├── Visualisation/
    │   ├── Create_Excel.py
    │   ├── output.xlsx
    │   ├── visunetcdf.py
    │ 
    ├── requirements.txt
  ```

## Description des Répertoires et Fichiers

### Dat_processing
Ce répertoire contient des scripts et dossiers pour le traitement des données liées à l'albédo, aux précipitations et à la température. Les scripts génèrent diverses figures et images basées sur les données traitées. Pour utiliser les fonctions il faut remplir dans le main le bon Path (les données sont dans Data ), ainsi que les dates de debut et de fin de l'evenement à observer. 

- **Data/**: Contient les .dat à utilisés pour la période.
- **albedo_figures/**: Contient des figures liées à l'albédopar jour (valeure journalière).
- **precip_albedo_figures/**: Contient des figures liées aux précipitations et de l'albédo par jour.
- **precip_temp_figures/**: Contient des figures liées aux précipitations et à la température par jour.
- **temperature_figures/**: Contient des figures liées à la température.
- **daily_albedo.csv**: Fichier CSV contenant les données d'albédo quotidien, utilisé pour créer un fichier Excel avec des images.
- **Read_Dat.py**: Script pour lire le fichier .dat de données entre la date de debut et la date de fin.
- **Read_DaybyDay.py**: Script pour lire le fichier .dat avec des figures pour chaque jour entre la date de début et de fin passée en paramètre.

### Melting_Layer
Ce répertoire contient des scripts pour détecter et tracer les couches de fusion en utilisant les données netcdf du micro rain radar.

- **MRR_ML_detection_auto.py**: Script pour la détection automatique des couches de fusion.

- **MRR_ML_detection.py**: Script pour la détection manuelle des couches de fusion.

- **MRR_ploting.py**: Script pour tracer les données les graphes des couches de fusion.

### Micro_Rain_Radar
Ce répertoire contient des scripts et des données liées au traitement du micro-radar de pluie.

- **.ipynb_checkpoints/**: Contient les checkpoints pour les notebooks Jupyter.

- **Data/**: Contient les fichiers de données brutes.

- **lib/**: Contient les fichiers de bibliothèque.
- **MK_processed/**: Répertoire pour stocker les images traitées.
  - **202309/**: Images traitées pour septembre 2023.
  - **202310/**: Images traitées pour octobre 2023.
  - **202311/**: Images traitées pour novembre 2023.
  - **PYRAMIDE/**: Images traitées pour le projet PYRAMIDE.
- **default_parameters.txt**: Contient les paramètres par défaut pour le traitement => ce sont les paramètres a modifiée avec la bonne PATH, ainsi que par exemple la height des images MRR.
- **requirements.txt**: Contient la liste des packages Python requis.
- **Run_MK.py**: Script principal pour exécuter le traitement MK.

### Visualisation
Ce répertoire contient la sortie de visualisation finale.

- **Create_Excel.py**: Script Permettant de Visualiser les données du MRR, mais aussi les données Dat si les dates correspondent. 

## Utilisation

1. **Installer les Dépendances**: Assurez-vous d'avoir tous les packages Python requis en exécutant :
    ```bash
    pip install -r requirements.txt
    ```

2. **Traitement des Données**: Utilisez les scripts dans les répertoires `Dat_processing` et `Melting_Layer` pour traiter les données et générer des figures.

3. **Exécuter le Traitement MK**: Exécutez le script `Run_MK.py` pour traiter les données du micro-radar de pluie et stocker les images dans le répertoire `MK_processed`.

4. **Générer la Visualisation**: Utilisez le script `Create_Excel.py` dans le répertoire `Visualisation` pour générer un fichier Excel avec des images.

## Contribuer
Si vous souhaitez contribuer à ce projet, veuillez suivre ces étapes :
1. Forkez le dépôt.
2. Créez une nouvelle branche (`git checkout -b feature-branch`).
3. Apportez vos modifications.
4. Commitez vos modifications (`git commit -m 'Ajouter une fonctionnalité'`).
5. Poussez sur la branche (`git push origin feature-branch`).
6. Ouvrez une Pull Request.

## Licence
Ce projet est sous licence MIT - voir le fichier LICENSE pour plus de détails.

## Remerciements
Merci à tous les contributeurs et collaborateurs qui ont rendu ce projet possible.(Geremy Panthou, Brice Boudevillain, Catherine Coulaud, Arno Reboud, Félix Delorme, Thomas Condom). 
