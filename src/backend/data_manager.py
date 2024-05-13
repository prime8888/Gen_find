import os
from Bio import Entrez, SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.Seq import Seq
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from urllib.error import HTTPError
from datetime import datetime

email = "daniil.kudriashov@etu.unistra.fr"

def fetch_nc_ids_for_organisms(organisms):
    Entrez.email = email
    retries = 5
    delay = 1  # start with 1 second delay
    updated_organisms = {}

    for organism, item in organisms.items():
        path, mod_date = item
        last_exception = None
        
        # Create a basic query for organism
        query = f"{organism}[Orgn] AND (\"NC_000000\"[Accession] : \"NC_999999\"[Accession])"

        # Add date filter if modification date is provided and not None
        if mod_date is not None:
            formatted_date = datetime.strptime(mod_date, '%Y-%m-%d %H:%M:%S')
            formatted_date = formatted_date.strftime('%Y/%m/%d')
            date_filter = f" AND (\"{formatted_date}\"[PDAT] : \"3000/01/01\"[PDAT])"  # Assuming future upper bound
            query += date_filter
        
        for attempt in range(retries):
            try:
                with Entrez.esearch(db="nucleotide", term=query, rettype="gb", retmode="text", retmax=9999999) as handle:
                    search_results = Entrez.read(handle)
                    if 'IdList' in search_results and search_results['IdList']:
                        print(f"Found {len(search_results['IdList'])} newer records for {organism}")
                        updated_organisms[organism] = {'ids': search_results['IdList'], 'path': path}
                        break  # Successfully fetched, break out of retry loop
                    else:
                        print(f"No newer records found for {organism}")
                        break
            except HTTPError as e:
                print(f"HTTP error fetching IDs for {organism}, attempt {attempt + 1}: {e}")
                last_exception = e
                if e.code == 429 or e.code == 503:
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                else:
                    break  # Break on other HTTP errors
            except Exception as exc:
                print(f"An error occurred fetching IDs for {organism}: {exc}")
                last_exception = exc
                break  # Exit the loop if an unexpected error occurs
        
        if last_exception and organism not in updated_organisms:
            print(f"Failed to fetch NCBI IDs for {organism} after {retries} attempts: {last_exception}")

    return updated_organisms

def fetch_genbank_records_by_organism(organism, data, rettype="gb"):
    """Fetch GenBank records for all IDs of a single organism with retries and exponential backoff.
    Skip the organism if all retries fail due to recoverable errors."""
    Entrez.email = email
    retries = 5
    delay = 1  # start with 1 second delay
    last_exception = None

    for attempt in range(retries):
        try:
            with Entrez.efetch(db="nucleotide", id=','.join(data['ids']), rettype=rettype, retmode="text") as handle:
                print(f"Fetched records for {organism}")
                records = list(SeqIO.parse(handle, "genbank"))
                return organism, {'records': records, 'path': data['path']}
        except HTTPError as e:
            last_exception = e
            if e.code == 429 or e.code == 503:  # Too many requests or service unavailable
                print(f"Rate limit hit, retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
                if delay > 60:  # Prevent excessively long wait times
                    break
            else:
                break  # Break on other HTTP errors as they might not be recoverable
        except Exception as exc:
            last_exception = exc
            print(f"An error occurred fetching records for {organism}: {exc}")
            break  # Exit the loop if an unexpected error occurs

    # If all retries are exhausted or a severe error occurs, log and skip this organism
    print(f"Failed to fetch records for {organism} after {retries} attempts: {last_exception}")
    return organism, {'records': [], 'path': data['path']}  # Return an empty record list for this organism

def fetch_genbank_record_by_id(organism, id, path, rettype="gb"):
    """Fetch GenBank record for a single ID with retries. Skip the record if all retries fail."""
    Entrez.email = email
    retries = 5
    delay = 1  # start with 1 second delay
    last_exception = None

    for attempt in range(retries):
        try:
            with Entrez.efetch(db="nucleotide", id=id, rettype=rettype, retmode="text") as handle:
                print(f"Fetched record for ID {id} in {organism}")
                record = list(SeqIO.parse(handle, "genbank"))
                print(f"Parsed record for ID {id} in {organism}")
                return organism, {'record': record, 'path': path}
        except HTTPError as e:
            if e.code == 429 or e.code == 503:  # Too many requests or service unavailable
                print(f"Rate limit hit, retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
                last_exception = e
            else:
                last_exception = e
                break  # Break on unrecoverable HTTP errors
        except Exception as e:
            last_exception = e
            time.sleep(delay)
            delay *= 2  # Exponential backoff on general errors
            if delay > 60:  # Prevent excessively long wait times
                break

    # If all retries are exhausted, log and skip this record
    print(f"Failed to fetch record for ID {id} in {organism} after {retries} attempts: {last_exception}")
    return organism, {'record': [], 'path': path}  # Return empty list for this ID's record

def fetch_genbank_records(organisms, max_workers=5, rettype="gb"):
    """Download GenBank records using multithreading, switching strategy based on number of organisms."""
    records = {org: {'records': [], 'path': data['path']} for org, data in organisms.items()}
    
    if len(organisms) > max_workers:
        # Use batch fetching by organism
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(fetch_genbank_records_by_organism, organism, data, rettype): organism for organism, data in organisms.items()}
            for future in as_completed(futures):
                organism, result = future.result()
                records[organism] = result
    else:
        # Use fetching by individual IDs
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for organism, data in organisms.items():
                for id in data['ids']:
                    futures.append(executor.submit(fetch_genbank_record_by_id, organism, id, data['path'], rettype))
            for future in as_completed(futures):
                organism, result = future.result()
                records[organism]['records'].extend(result['record'])

    return records

def extract_cds_and_related_features(records, selected_regions):
    """Extracts CDS and, if selected, their exons and introns."""
    for organism, data in records.items():
        path = data['path']
        organism_modified = organism.replace(" ", "_")
        for record in data['records']:
            cds_data = {}
            introns_data = []
            if record:
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
            if record:
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

    sequence_length = len(record.seq)

    # Check if the location is a CompoundLocation to handle joins and complements
    if isinstance(feature.location, CompoundLocation):
        previous_end = None
        # Iterate through parts of the CompoundLocation
        for i, part in enumerate(feature.location.parts):
            # Validate bounds before processing
            if not validate_bounds(part.start, part.end, sequence_length):
                continue  # Skip this part if it is out of bounds
            
            # Extract the sequence for each part
            part_seq = record.seq[part.start:part.end]
            # If the part is on the reverse strand, reverse complement it
            if 'complement' in feature.location.operator:
                part_seq = part_seq.reverse_complement()
            sequence += part_seq
            exons.append({'sequence': part_seq, 'location': f"{part.start}..{part.end}"})

            # Calculate introns
            if previous_end is not None and part.start > previous_end + 1:
                if validate_bounds(previous_end + 1, part.start - 1, sequence_length):  # Validate intron bounds
                    intron_seq = record.seq[previous_end:part.start]
                    if 'complement' in feature.location.operator:
                        intron_seq = intron_seq.reverse_complement()
                    introns.append({'sequence': intron_seq, 'location': f"{previous_end+1}..{part.start-1}", 'start_exon': i, 'end_exon': i + 1})
            previous_end = part.end

        # If the overall operator is complement, reverse complement the entire sequence
        if 'complement' in feature.location.operator:
            sequence = sequence.reverse_complement()
    else:
        # For a simple location, validate bounds before extracting the sequence
        start = int(feature.location.start)
        end = int(feature.location.end)
        if validate_bounds(start, end, sequence_length):
            sequence = feature.location.extract(record.seq)
        else:
            # Handle the case where the location is out of bounds
            sequence = Seq('')  # Return an empty sequence if out of bounds

    return sequence, exons, introns

def validate_bounds(start, end, sequence_length):
    """Validate bounds of a genomic region."""
    return start < end and 0 <= start < sequence_length and 0 < end <= sequence_length