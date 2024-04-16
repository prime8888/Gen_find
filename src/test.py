from Bio import Entrez
from Bio import SeqIO

# Toujours fournir votre adresse e-mail réelle
Entrez.email = "daniil.kudriashov@etu.unistra.fr"

# Recherche des séquences génomiques pour Homo sapiens
# search_query = "Homo sapiens[orgn]"
# database = "nucleotide"
# max_records = 10  # Le nombre maximum d'enregistrements que vous souhaitez récupérer

# search_handle = Entrez.esearch(db=database, term=search_query, retmax=max_records)
# search_results = Entrez.read(search_handle)
# search_handle.close()

# # ID list contient les identifiants de séquence pour Homo sapiens
# id_list = search_results['IdList']


# # Récupération des détails des séquences à l'aide des identifiants
# fetch_handle = Entrez.efetch(db=database, id=id_list, rettype="gb", retmode="text")
# data = fetch_handle.read()
# fetch_handle.close()

# # Sauvegarde des données dans un fichier
# with open("homo_sapiens_sequences.gb", "w") as file:
#     file.write(data)

# Pour chaque enregistrement dans le fichier GenBank
for record in SeqIO.parse("homo_sapiens_sequences.gb", "genbank"):
    print(f"Analyse de l'enregistrement: {record.id}")
    for feature in record.features:
        if feature.type == "CDS":
            # Extraire la position de début et de fin
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            cds_sequence = feature.extract(record.seq)
            
            # Afficher les informations
            print(f"Identifiant CDS: {feature.qualifiers.get('protein_id', [''])[0]}")
            print(f"Position de début: {start}, Position de fin: {end}, Brin: {strand}")
            print(f"Séquence CDS: {cds_sequence}")
            
            # Si la fonctionnalité contient une séquence de protéines traduite
            if 'translation' in feature.qualifiers:
                protein_sequence = feature.qualifiers['translation'][0]
                print(f"Séquence de protéine: {protein_sequence}")
                
            # Sauvegarde des informations
            with open(f"{record.id}_cds.txt", "w") as output_file:
                output_file.write(f"> {record.id}_cds\n{cds_sequence}\n")
            with open(f"{record.id}_protein.txt", "w") as output_file:
                output_file.write(f"> {record.id}_protein\n{protein_sequence}\n")