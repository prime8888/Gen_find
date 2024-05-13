from backend import *

def main():
    # Execution
    # fetch_and_update_overview()
    # get organisms from path
    organisms = get_organisms_from_path("C:/Educ/Analyse_texte/new_proj/Results/Eukaryota/Animals/Mammals")
    print(organisms)
    organisms = fetch_nc_ids_for_organisms(organisms)
    records = fetch_genbank_records(organisms)
    print("Records fetched. Extracting...")
    process_genomic_data(records, ["CDS", "intron", "mRNA"])
    # extract_and_save_regions(records)



if __name__ == "__main__":
    main()