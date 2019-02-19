Cration d'un package conda
==========================

1. Créer un recipe et d'un fichier  :

Le recipe est composé de 3 scripts :

	* meta.yaml contient les metadata du recipe, seulement le nom et la version du package sont obligatoires
	* build.sh contient les commandes pour construire le package sur macOS et Linux
	* bld.bat contient les commandes pour construire le package sur Windows

2. Céer un fichier setup.py

setup.py est exécuté par les scripts builssh et bld.bat lors de la construction du package

3. Construire le package et l'upload sur le cloud anaconda

	conda build chemin_du_recipe (création d'un fichier tar.bz2)
	anaconda upload chemin_du_fichier_tar.bz2 (créer un compte sur le cloud avant)

4. Créer un environnement conda

	conda create -n nom_env
	conda activate
	conda activate nom_env

5. Installer le package

	conda config --add channels nom_channel
	conda install nom_package

6. Documentation

* Création du récipé : [https://conda.io/projects/conda-build/en/latest/source/recipe.html](https://conda.io/projects/conda-build/en/latest/source/recipe.html)
* Création du meta.yaml : [https://conda.io/projects/conda-build/en/latest/source/resources/define-metadata.html](https://conda.io/projects/conda-build/en/latest/source/resources/define-metadata.html)
* Création du setup.py : 
	* [https://docs.python.org/2/distutils/setupscript.html](https://docs.python.org/2/distutils/setupscript.html)
	* [https://packaging.python.org/tutorials/packaging-projects/](https://packaging.python.org/tutorials/packaging-projects/)
	* [https://deusyss.developpez.com/tutoriels/Python/packaging_pypi/](https://deusyss.developpez.com/tutoriels/Python/packaging_pypi/)
* Construction du package conda : [https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages/#cloud-uploading-conda-packages](https://docs.anaconda.com/anaconda-cloud/user-guide/tasks/work-with-packages/#cloud-uploading-conda-packages)
* Cloud anaconda : [https://anaconda.org/](https://anaconda.org/)








