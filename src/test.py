from Bio import Entrez

# Toujours fournir votre adresse e-mail réelle
Entrez.email = "daniil.kudriashov@etu.unistra.fr"

# Recherche des séquences génomiques pour Homo sapiens
search_query = "Homo sapiens[orgn]"
database = "nucleotide"
max_records = 10  # Le nombre maximum d'enregistrements que vous souhaitez récupérer

search_handle = Entrez.esearch(db=database, term=search_query, retmax=max_records)
search_results = Entrez.read(search_handle)
search_handle.close()

# ID list contient les identifiants de séquence pour Homo sapiens
id_list = search_results['IdList']


# Récupération des détails des séquences à l'aide des identifiants
fetch_handle = Entrez.efetch(db=database, id=id_list, rettype="gb", retmode="text")
data = fetch_handle.read()
fetch_handle.close()

# Sauvegarde des données dans un fichier
with open("homo_sapiens_sequences.gb", "w") as file:
    file.write(data)