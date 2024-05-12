from backend import *

def main():
    # Execution
    # fetch_and_update_overview()
    # get current path
    current_path = os.path.dirname(os.path.realpath(__file__))
    # get organisms from path
    organisms = get_organisms_from_path("C:/Educ/Analyse_texte/new_proj/Results/Eukaryota/Animals/Mammals/Homo sapiens")
    print(organisms)
    organisms = fetch_nc_ids_for_organisms(organisms)
    print(organisms)
    records = fetch_genbank_records(organisms)
    # print(records)
    print("Records fetched. Extracting...")
    extract_and_save_regions(records, ["CDS"])
    # extract_and_save_regions(records)



if __name__ == "__main__":
    main()