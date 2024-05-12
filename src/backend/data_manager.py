import os
from Bio import Entrez, SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.Seq import Seq

email = "daniil.kudriashov@etu.unistra.fr"

def fetch_nc_ids_for_organisms(organisms):
    """Fetch NCBI identifiers for specified organisms."""
    Entrez.email = email
    for organism, path in organisms.items():
        query = f"{organism}[Orgn] AND refseq[filter] AND (\"NC_000000\"[Accession] : \"NC_999999\"[Accession])"
        handle = Entrez.esearch(db="nucleotide", term=query, rettype="gb", retmode="text", retmax=9999999)
        search_results = Entrez.read(handle)
        handle.close()
        organisms[organism] = {'ids': search_results['IdList'], 'path': path}
    return organisms

def fetch_genbank_records(organisms):
    """ Télécharge les enregistrements GenBank pour une liste d'identifiants NC. """
    Entrez.email = email
    records = {}
    
    for organism, data in organisms.items():
        records[organism] = {'records': [], 'path': data['path']}
        for nc_id in data['ids']:
            handle = Entrez.efetch(db="nucleotide", id=nc_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            records[organism]['records'].append(record)
    return records

def extract_and_save_regions(records, selected_regions):
    """Extracts specified functional regions and saves them to files."""
    for organism, data in records.items():
        path = data['path']
        for record in data['records']:
            for feature in record.features:
                if feature.type in selected_regions:
                    # qualifiers = feature.qualifiers
                    # organism = qualifiers.get("organism", ["Unknown"])[0]
                    organism = organism.replace(" ", "_")
                    feature_type = feature.type
                    start = int(feature.location.start)
                    end = int(feature.location.end)

                    # Handle complex locations
                    sequence = handle_complex_location(feature, record)

                    if validate_bounds(start, end, len(record.seq)):
                        filename = f"{feature_type}_{organism}_{record.id}.txt"
                        file_path = os.path.join(path, filename)
                        
                        with open(file_path, 'w') as file:
                            file.write(f">{feature_type} {organism} NC_{record.id}: {start}..{end}\n")
                            file.write(str(sequence) + "\n")
                        print(f"Data saved to {file_path}")
                    else:
                        print(f"Invalid bounds for {feature_type} in NC_{record.id}: {start}-{end}")

def handle_complex_location(feature, record):
    """Extract sequence considering complex location types like join and complement."""
    sequence = Seq('')  # Initialize an empty Seq object

    # Check if the location is a CompoundLocation to handle joins and complements
    if isinstance(feature.location, CompoundLocation):
        # Iterate through parts of the CompoundLocation
        for part in feature.location.parts:
            # Extract the sequence for each part
            part_seq = record.seq[part.start:part.end]
            # If the part is on the reverse strand, reverse complement it
            if part.strand == -1:
                part_seq = part_seq.reverse_complement()
            sequence += part_seq

        # If the overall operator is complement, reverse complement the entire sequence
        if 'complement' in feature.location.operator:
            sequence = sequence.reverse_complement()
    else:
        # If it's a simple location, just extract the sequence
        sequence = feature.location.extract(record.seq)

    return sequence

def validate_bounds(start, end, sequence_length):
    """ Valide les bornes de la région fonctionnelle. """
    return start < end and 0 <= start < sequence_length and 0 < end <= sequence_length