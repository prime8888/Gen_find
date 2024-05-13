from backend import *
import os

def main():
    # Execution
    # fetch_and_update_overview()
    # get organisms from path
    organisms = get_organisms_from_path("C:/Educ/Analyse_texte/new_proj/Results/Eukaryota/Animals/Mammals/Hippopotamus amphibius")
    print(organisms)
    max_workers = os.cpu_count() + 2
    organisms = fetch_nc_ids_for_organisms(organisms)
    records = fetch_genbank_records(organisms, rettype="gbwithparts", max_workers=max_workers)
    print("Records fetched. Extracting...")
    process_genomic_data(records, ["CDS", "intron", "mRNA"])
    # extract_and_save_regions(records)



if __name__ == "__main__":
    main()