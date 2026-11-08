"""
Step 1 - Retrieve human alternatively spliced proteins from UniProt and
keep entries with >= 2 PDB structures deposited after the AlphaFold2
training cut-off (30 April 2018).

Code excerpt as published in Chapter 4.1 of the thesis. Helper functions
not printed in the thesis (e.g. `process_protein`) are marked TODO.
"""
import json
import multiprocessing as mp
import os
import time
from collections import defaultdict
from datetime import datetime

import requests
from tqdm import tqdm

UNIPROT_QUERY = (
    "https://rest.uniprot.org/uniprotkb/stream?format=json&query="
    "var_seq%20AND%20%28model_organism%3A9606%29%20AND%20"
    "%28proteins_with%3A1%29%20AND%20%28proteins_with%3A5%29"
)


def get_protein_data():
    if os.path.exists("./uniprot_data.json"):
        with open("./uniprot_data.json", "r") as fh:
            return json.load(fh)
    response = requests.get(UNIPROT_QUERY)
    if response.status_code == 200:
        data = response.json()
        with open("./uniprot_data.json", "w") as fh:
            json.dump(data, fh, indent=4)
        return data
    return None


def filter_relevant_data(protein_data):
    filtered_data = defaultdict(list)
    for protein in protein_data:
        accession_id = protein["primaryAccession"]
        for comment in protein.get("comments", []):
            if "isoforms" in comment:
                filtered_data[accession_id].extend(comment["isoforms"])
    return filtered_data


def get_pdb_deposition_date(pdb_id, session):
    cache_file = f"./pdb_cache/{pdb_id}.json"
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            return datetime.fromisoformat(json.load(f)["date"])

    for _ in range(3):  # retry mechanism
        try:
            api_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            response = session.get(api_url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                deposition_date = data.get("rcsb_accession_info", {}).get("deposit_date")
                if deposition_date:
                    date_obj = datetime.fromisoformat(deposition_date.split("T")[0])
                    os.makedirs("./pdb_cache", exist_ok=True)
                    with open(cache_file, "w") as f:
                        json.dump({"date": date_obj.isoformat()}, f)
                    return date_obj
            time.sleep(1)
        except Exception as e:
            print(f"Error processing PDB ID {pdb_id}: {e}")
            time.sleep(1)
    return None


def process_protein(protein):
    """TODO: original implementation not printed in the thesis.
    Returns (accession_id, [PDB IDs deposited after 2018-04-30])."""
    raise NotImplementedError


def pdb_data(protein_data):
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(process_protein, protein_data),
                            total=len(protein_data), desc="Processing proteins"))
    return {acc_id: structs for acc_id, structs in results if structs}


def main():
    protein_data = get_protein_data()
    if protein_data and "results" in protein_data:
        filtered_protein_data = filter_relevant_data(protein_data["results"])
        pdb_id_data = pdb_data(protein_data["results"])

        relevant_data = {
            acc: {"isoforms": isoforms, "pdb_structures": pdb_id_data[acc]}
            for acc, isoforms in filtered_protein_data.items()
            if acc in pdb_id_data and len(pdb_id_data[acc]) >= 2
        }
        if relevant_data:
            with open("./filtered_protein_with_pdb.json", "w") as fh:
                json.dump(relevant_data, fh, indent=4)


if __name__ == "__main__":
    main()
