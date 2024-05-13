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
    """Returns a dictionary of organism names with their corresponding paths and most recent modification dates by traversing directories recursively."""
    organisms = {}
    
    def traverse_directories(current_path):
        entries = os.listdir(current_path)
        most_recent_date = None
        formatted_date = None
        
        for entry in entries:
            path = os.path.join(current_path, entry)
            if os.path.isdir(path):
                # Continue to traverse the directory
                traverse_directories(path)
            else:
                # Get the modification time for each file
                modification_time = os.path.getmtime(path)
                if most_recent_date is None or modification_time > most_recent_date:
                    most_recent_date = modification_time
        
        if most_recent_date is not None:
            # Format the most recent modification date
            formatted_date = datetime.fromtimestamp(most_recent_date).strftime('%Y-%m-%d %H:%M:%S')
        organism_name = os.path.basename(current_path)
            
        organisms[organism_name] = (current_path, formatted_date)

    traverse_directories(selected_path)
    return organisms

def get_organisms_from_path_list(selected_paths):
    """Returns a dictionary of organism names with their corresponding paths and most recent modification dates from a list of paths."""
    organisms = {}
    for path in selected_paths:
        orgs = get_organisms_from_path(path)
        organisms.update(orgs)
    return organisms