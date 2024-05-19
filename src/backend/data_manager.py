import os
from Bio import Entrez, SeqIO
from Bio.SeqFeature import CompoundLocation
from Bio.Seq import Seq
import time
import random
from urllib.error import HTTPError
from datetime import datetime
from ui.logger import *

email = "daniil.kudriashov@etu.unistra.fr"

def fetch_nc_ids_for_organisms(organisms):
    Entrez.email = email
    retries = 5
    delay = 1  # start with 1 second delay
    updated_organisms = {}
    logger = logging.getLogger()

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
                        logger.info(f"Found {len(search_results['IdList'])} newer records for {organism}")
                        updated_organisms[organism] = {'ids': search_results['IdList'], 'path': path}
                        break  # Successfully fetched, break out of retry loop
                    else:
                        logger.info(f"No newer records found for {organism}")
                        break
            except HTTPError as e:
                logger.warning(f"HTTP error fetching IDs for {organism}, attempt {attempt + 1}: {e}")
                last_exception = e
                if e.code == 429 or e.code == 503:
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                else:
                    break  # Break on other HTTP errors
            except Exception as exc:
                logger.warning(f"An error occurred fetching IDs for {organism}: {exc}")
                last_exception = exc
                break  # Exit the loop if an unexpected error occurs
        
        if last_exception and organism not in updated_organisms:
            logger.warning(f"Failed to fetch NCBI IDs for {organism} after {retries} attempts: {last_exception}")

    return updated_organisms


def fetch_records_and_process(organism, id, path, selected_regions, rettype="gbwithparts"):
    """Fetch GenBank record for a single ID with retries. Skip the record if all retries fail."""
    Entrez.email = email
    retries = 5
    delay = 1  # start with 1 second delay
    last_exception = None
    logger = logging.getLogger()

    time.sleep(random.uniform(0, 1))
    for attempt in range(retries):
        try:
            if(attempt > 0):
                logger.warning(f"Attempt {attempt + 1} to fetch record for ID {id} in {organism}")
            with Entrez.efetch(db="nucleotide", id=id, rettype=rettype, retmode="text") as handle:
                record = list(SeqIO.parse(handle, "genbank"))
                data = {organism: {'records': record, 'path': path}}
                process_genomic_data(data, selected_regions)
                logger.info(f"Processed record for ID {id} in {organism}")
                return
        except HTTPError as e:
            if e.code == 429 or e.code == 503:  # Too many requests or service unavailable
                logger.warning(f"Rate limit hit, retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
                last_exception = e
            else:
                last_exception = e
                logger.error(f"Exception for ID {id} in {organism}: {e}")
                break  # Break on unrecoverable HTTP errors
        except Exception as e:
            last_exception = e
            time.sleep(delay)
            delay *= 2  # Exponential backoff on general errors
            if delay > 60:  # Prevent excessively long wait times
                break
            

    # If all retries are exhausted, log and skip this record
    logger.warning(f"Failed to fetch record for ID {id} in {organism} after {retries} attempts: {last_exception}")
    return


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
                    feature_loc = ''
                    if isinstance(feature.location, CompoundLocation) and "join" in str(feature.location) and feature.location.strand == -1:
                        feature_loc = "join("
                        for part in reversed(feature.location.parts):
                            feature_loc += str(part) + ","
                        feature_loc = feature_loc[:-1] + ")"
                    else:
                        feature_loc = str(feature.location)
                    if feature.location.strand == -1:
                        feature_loc = "complement(" + feature_loc + ")"
                    if feature.type == "CDS" and feature_loc not in cds_data:
                        sequence, exons, introns = handle_complex_location(feature, record)
                        cds_loc = feature_loc.replace("(+)", "").replace("(-)", "").replace(":", "..").replace("{", "(").replace("}", ")").replace("[", "").replace("]", "").replace(" ", "")
                        cds_info = {
                            'header': f"CDS {organism} {record.name}: {cds_loc}",
                            'sequence': sequence,
                            'location': cds_loc,
                            'exons': exons,
                            'introns': introns
                        }
                        cds_data[feature_loc] = cds_info

                # Write CDS data and exon details if selected
                if cds_data:
                    filename = f"CDS_{organism_modified}_{record.name}.txt"
                    file_path = os.path.join(path, filename)
                    with open(file_path, 'w') as file:
                        for cds in cds_data.values():
                            if len(cds['sequence']) > 0:
                                file.write(f"{cds['header']}\n")
                                file.write(str(cds['sequence']) + "\n")
                            for i, exon in enumerate(cds['exons'], start=1):
                                if len(exon['sequence']) > 0:
                                    file.write(f"CDS_{organism_modified} {record.name}: {cds['location']} Exon {i}\n")
                                    file.write(str(exon['sequence']) + "\n")
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
                            if len(intron['sequence']) > 0:
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
                    feature_loc = ''
                    if isinstance(feature.location, CompoundLocation) and "join" in str(feature.location) and feature.location.strand == -1:
                        feature_loc = "join("
                        for part in reversed(feature.location.parts):
                            feature_loc += str(part) + ","
                        feature_loc = feature_loc[:-1] + ")"
                    else:
                        feature_loc = str(feature.location)
                    if feature.location.strand == -1:
                        feature_loc = "complement(" + feature_loc + ")"
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
                # Write collected data for each type of region to corresponding files
                for feature_type, items in all_region_data.items():
                    filename = f"{feature_type}_{organism_modified}_{record.name}.txt"
                    file_path = os.path.join(path, filename)
                    if len(items) > 0:
                        with open(file_path, 'w') as file:
                            for region in items.values():
                                if len(region['sequence']) > 0:
                                    file.write(f"{region['header']}\n")
                                    file.write(str(region['sequence']) + "\n")

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
        previous_start = None
        # Iterate through parts of the CompoundLocation
        for i, part in enumerate(feature.location.parts):
            # Validate bounds before processing
            if not validate_bounds(part.start, part.end, sequence_length):
                continue  # Skip this part if it is out of bounds
            
            # Extract the sequence for each part
            part_seq = record.seq[part.start:part.end]
            # If the part is on the reverse strand, reverse complement it
            if feature.location.strand == -1:
                part_seq = part_seq.reverse_complement()
            sequence += part_seq
            exons.append({'sequence': part_seq, 'location': f"{part.start}..{part.end}"})

            # Calculate introns
            if feature.location.strand == -1:
                if previous_start is not None:
                    if validate_bounds(part.end + 1, previous_start - 1, sequence_length):  # Validate intron bounds
                        intron_seq = record.seq[part.end:previous_start]
                        intron_seq = intron_seq.reverse_complement()
                        introns.append({'sequence': intron_seq, 'location': f"{part.end+1}..{previous_start-1}", 'start_exon': i, 'end_exon': i + 1})
                previous_start = part.start
            else:
                if previous_end is not None:
                    if validate_bounds(previous_end + 1, part.start - 1, sequence_length):  # Validate intron bounds
                        intron_seq = record.seq[previous_end:part.start]
                        introns.append({'sequence': intron_seq, 'location': f"{previous_end+1}..{part.start-1}", 'start_exon': i, 'end_exon': i + 1})
                previous_end = part.end

    else:
        # For a simple location, validate bounds before extracting the sequence
        start = int(feature.location.start)
        end = int(feature.location.end)
        if validate_bounds(start, end, sequence_length):
            sequence = record.seq[feature.location.start:feature.location.end]
            if feature.location.strand == -1:
                sequence = sequence.reverse_complement()
        else:
            # Handle the case where the location is out of bounds
            sequence = Seq('')  # Return an empty sequence if out of bounds
    return sequence, exons, introns

def validate_bounds(start, end, sequence_length):
    """Validate bounds of a genomic region."""
    return start < end and 0 <= start < sequence_length and 0 < end <= sequence_length