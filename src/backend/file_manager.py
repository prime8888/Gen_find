import requests
import pandas as pd
import os
from datetime import datetime

url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"
local_path = "local_overview.txt"

def fetch_and_update_overview():
    """Fetch overview file from NCBI FTP and check if updates are needed based on the 'last-modified' header."""
    if os.path.exists(local_path):
        local_modified_time = datetime.fromtimestamp(os.path.getmtime(local_path))
    else:
        local_modified_time = datetime.min

    response = requests.head(url)
    last_modified = response.headers.get('Last-Modified')
    print(f"Last modified time: {last_modified}")
    last_modified_time = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S GMT')
    print(f"Last modified time: {last_modified_time}")

    if last_modified_time > local_modified_time or not os.path.exists(local_path):
        print("Newer file found or local file does not exist, downloading update...")
        response = requests.get(url)
        with open(local_path, 'wb') as file:
            file.write(response.content)
        create_directories_from_overview()
    else:
        print("Local data is up-to-date.")

def create_directories_from_overview():
    """Creates directories based on the organism information in the overview file."""
    df = pd.read_csv(local_path, sep='\t')
    base_path = "./Results"

    for index, row in df.iterrows():
        directory_path = os.path.join(base_path, row['Kingdom'], row['Group'], row['SubGroup'], row['#Organism/Name'].replace("[", "").replace(":", "").replace("?", "").replace("/", "").replace("]", ""))
        os.makedirs(directory_path, exist_ok=True)

def get_organisms_from_path(selected_path):
    """Returns a dictionary of organism names with their corresponding paths by traversing directories recursively."""
    organisms = {}
    # Fonction récursive pour parcourir les sous-répertoires
    def traverse_directories(current_path):
        # Liste les entrées dans le répertoire courant
        entries = os.listdir(current_path)
        if any(os.path.isdir(os.path.join(current_path, entry)) for entry in entries):
            # Parcours des sous-répertoires s'il y en a
            for entry in entries:
                path = os.path.join(current_path, entry)
                if os.path.isdir(path):
                    traverse_directories(path)
        else:
            # Ajoute le dernier niveau de répertoire (organisme) au dictionnaire avec le chemin complet
            organism_name = os.path.basename(current_path)
            organisms[organism_name] = current_path

    # Commence la traversée depuis le chemin sélectionné
    traverse_directories(selected_path)
    return organisms
