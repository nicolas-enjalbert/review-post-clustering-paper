#!/bin/bash

# Dossier contenant les scripts R
SCRIPT_DIR="scripts"

# Vérifier si le dossier existe
if [ ! -d "$SCRIPT_DIR" ]; then
    echo "Le dossier $SCRIPT_DIR n'existe pas."
    exit 1
fi

# Boucle sur tous les fichiers .R du dossier
for script in "$SCRIPT_DIR"/*.R; do
    echo "Exécution de $script ..."
    R --vanilla < "$script" > "${script%.R}.log" 2>&1 &
done

echo "Tous les scripts sont lancés en arrière-plan."

