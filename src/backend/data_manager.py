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
        handle = Entrez.efetch(db="nucleotide", id=','.join(data['ids']), rettype="gbwithparts", retmode="text")
        print("Fetched records for", organism, "starting to parse...")
        records[organism]['records'] = list(SeqIO.parse(handle, "genbank"))
        handle.close()
    return records

def extract_cds_and_related_features(records, selected_regions):
    """Extracts CDS and, if selected, their exons and introns."""
    for organism, data in records.items():
        path = data['path']
        organism_modified = organism.replace(" ", "_")
        for record in data['records']:
            cds_data = {}
            introns_data = []
            for feature in record.features:
                feature_loc = str(feature.location)
                if feature.type == "CDS" and feature_loc not in cds_data:
                    sequence, exons, introns = handle_complex_location(feature, record)
                    cds_loc = feature_loc.replace("(+)", "").replace("(-)", "").replace(":", "..").replace("{", "(").replace("}", ")").replace("[", "").replace("]", "") 
                    cds_info = {
                        'header': f"CDS {organism} {record.name}: {cds_loc}",
                        'sequence': sequence,
                        'location': cds_loc,
                        'exons': exons,
                        'introns': introns
                    }
                    cds_data[feature_loc] = cds_info

            # Write CDS data and exon details if selected
            # print(cds_data)
            if cds_data:
                filename = f"CDS_{organism_modified}_{record.name}.txt"
                file_path = os.path.join(path, filename)
                with open(file_path, 'w') as file:
                    for cds in cds_data.values():
                        file.write(f"{cds['header']}\n")
                        file.write(str(cds['sequence']) + "\n")
                        for i, exon in enumerate(cds['exons'], start=1):
                            file.write(f"CDS_{organism_modified} {record.name}: {cds['location']} Exon {i}\n")
                            file.write(str(exon['sequence']) + "\n")
                        file.write("\n")
                        if 'intron' in selected_regions and cds['introns']:
                            for i, intron in enumerate(cds['introns'], start=1):
                                introns_data.append({
                                    'header': f"intron {organism} {record.name}: {cds['location']} Intron {i}",
                                    'sequence': intron['sequence']
                                })

            # Write intron data if introns are selected
            if 'intron' in selected_regions and introns_data:
                intron_filename = f"intron_{organism_modified}_{record.name}.txt"
                intron_file_path = os.path.join(path, intron_filename)
                with open(intron_file_path, 'w') as file:
                    for intron in introns_data:
                        file.write(f"{intron['header']}\n")
                        file.write(str(intron['sequence']) + "\n")

def extract_other_functional_regions(records, selected_regions):
    """Extracts and saves all specified functional regions (except CDS) to files."""
    for organism, data in records.items():
        path = data['path']
        organism_modified = organism.replace(" ", "_")
        for record in data['records']:
            all_region_data = {}  # Dictionary to store all regions data before writing
            for feature in record.features:
                feature_loc = str(feature.location)
                if feature.type in selected_regions and feature.type != "CDS" and (feature.type not in all_region_data or feature_loc not in all_region_data[feature.type]):
                    sequence, _, _ = handle_complex_location(feature, record)
                    feature_type = feature.type
                    
                    
                    # Collect data for each type of region
                    if feature_type not in all_region_data:
                        all_region_data[feature_type] = {}

                    
                    region_loc = feature_loc.replace("(+)", "").replace("(-)", "").replace(":", "..").replace("{", "(").replace("}", ")").replace("[", "").replace("]", "")

                    region_info = {
                        'header': f"{feature_type} {organism} {record.name}: {region_loc}",
                        'sequence': sequence
                    }
                    all_region_data[feature_type][feature_loc] = region_info
            # print(all_region_data)
            # Write collected data for each type of region to corresponding files
            for feature_type, items in all_region_data.items():
                filename = f"{feature_type}_{organism_modified}_{record.name}.txt"
                file_path = os.path.join(path, filename)
                if len(items) > 0:
                    with open(file_path, 'w') as file:
                        for region in items.values():
                                file.write(f"{region['header']}\n")
                                file.write(str(region['sequence']) + "\n\n")

def process_genomic_data(records, selected_regions):
    if "CDS" in selected_regions or "exon" in selected_regions or "intron" in selected_regions:
        extract_cds_and_related_features(records, selected_regions)
    extract_other_functional_regions(records, selected_regions)

def handle_complex_location(feature, record):
    """Extract sequence considering complex location types like join and complement."""
    sequence = Seq('')  # Initialize an empty Seq object
    exons = []
    introns = []

    # Check if the location is a CompoundLocation to handle joins and complements
    if isinstance(feature.location, CompoundLocation):
        previous_end = None
        # Iterate through parts of the CompoundLocation
        for i, part in enumerate(feature.location.parts):
            # Extract the sequence for each part
            part_seq = record.seq[part.start:part.end]
            # If the part is on the reverse strand, reverse complement it
            if 'complement' in feature.location.operator:
                part_seq = part_seq.reverse_complement()
            sequence += part_seq
            exons.append({'sequence': part_seq, 'location': f"{part.start}..{part.end}"})

            # Calculate introns
            if previous_end is not None and part.start > previous_end + 1:
                intron_seq = record.seq[previous_end:part.start]
                if 'complement' in feature.location.operator:
                    intron_seq = intron_seq.reverse_complement()
                introns.append({'sequence': intron_seq, 'location': f"{previous_end+1}..{part.start-1}", 'start_exon': i, 'end_exon': i + 1})
            previous_end = part.end

        # If the overall operator is complement, reverse complement the entire sequence
        if 'complement' in feature.location.operator:
            sequence = sequence.reverse_complement()
    else:
        # If it's a simple location, just extract the sequence
        sequence = feature.location.extract(record.seq)

    return sequence, exons, introns

def validate_bounds(start, end, sequence_length):
    """ Valide les bornes de la région fonctionnelle. """
    return start < end and 0 <= start < sequence_length and 0 < end <= sequence_length