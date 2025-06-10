import os
import random
import logging
import requests
from io import StringIO
from Bio import SeqIO

from app.helpers.Duru.core.sequence_processing import create_fasta_from_sequences, SequenceProcessingError

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

class BOLDAPIError(Exception):
    pass

class SequenceData:
    def __init__(self, processid, organism_name, marker_code, sequence, query_name):
        self.processid = processid
        self.organism_name = organism_name
        self.marker_code = marker_code
        self.sequence = sequence
        self.query_name = query_name

    def __repr__(self):
        return (f"SequenceData(processid='{self.processid}', "
                f"organism='{self.organism_name}', "
                f"marker='{self.marker_code}', "
                f"query_name='{self.query_name}', "
                f"sequence='{self.sequence[:30]}...')")  # Show only first 30 chars of seq

BOLD_API_URL = "http://www.boldsystems.org/index.php/API_Public"

MARKER_ALIASES = {
    "COI-5P": {"COI-5P", "CO1", "COI", "COI5P", "COI-5’"},
    "COI-3P": {"COI-3P", "COI3P", "COI-3’"},
    "COI": {"COI", "CO1"},  # Broad COI to catch variants not specifically 5P/3P
    "16S": {"16S", "16S rRNA", "16SrRNA"},
    "18S": {"18S", "18S rRNA", "18SrRNA"},
    "12S": {"12S", "12S rRNA", "12SrRNA"},
    "ITS": {"ITS", "ITS1", "ITS2", "INTERNAL TRANSCRIBED SPACER"},
    "rbcL": {"RBCL", "rbcL", "RIBP"},
    "matK": {"MATK", "matK"},
    "CYTB": {"CYTB", "CYTOCHROME B", "CYT B"},
    "ND2": {"ND2", "NADH DEHYDROGENASE SUBUNIT 2"},
    "ND4": {"ND4"},
    "trnH-psbA": {"TRNH-PSBA", "TRNHPSBA", "TRN-H PSBA"},
    "COII": {"COII", "COX2"},
    "COIII": {"COIII", "COX3"},
    "ATP6": {"ATP6"},
    "ATP8": {"ATP8"},
}

def get_canonical_marker(marker: str) -> str | None:
    normalized = marker.strip().upper()
    for canonical, aliases in MARKER_ALIASES.items():
        if normalized in {a.upper() for a in aliases}:
            return canonical
    return None

def is_matching_marker(user_marker: str, found_marker: str) -> bool:
    user_canonical = get_canonical_marker(user_marker)
    found_canonical = get_canonical_marker(found_marker)
    return user_canonical is not None and user_canonical == found_canonical

def fetch_sequences_by_organism_and_marker(organism: str, marker_code: str, limit: int = 10) -> list[SequenceData]:
    api_marker_code = marker_code.replace("-", "").upper()

    params = {
        'taxon': organism,
        'marker': api_marker_code,
    }

    url = f"{BOLD_API_URL}/sequence"

    logging.info(f"Requesting URL: {url} with params: {params}")

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        logging.info(f"API raw response snippet (first 500 chars):\n{response.text[:500]}")

        if not response.text.strip():
            logging.info("API returned empty response text.")
            return []
    except requests.exceptions.RequestException as e:
        raise BOLDAPIError(f"Failed to connect to BOLD Systems API: {e}")

    if response.text.strip().lower().startswith("<!doctype html>") or \
       response.text.strip().lower().startswith("<html"):
        raise BOLDAPIError(f"BOLD Systems API returned an HTML page (likely an error or 'no results' page), not FASTA. "
                           f"Partial HTML: {response.text[:200]}...")

    sequences = []
    try:
        for record in SeqIO.parse(StringIO(response.text), "fasta"):
            header_parts = record.id.split('|')
            processid = header_parts[0] if len(header_parts) > 0 else "N/A"
            organism_name_found = header_parts[1] if len(header_parts) > 1 else "N/A"

            marker_code_found = "N/A"
            if len(header_parts) > 2:
                marker_code_found = header_parts[2]
            elif record.description and "marker=" in record.description:
                desc_parts = record.description.split(' ')
                for part in desc_parts:
                    if part.startswith("marker="):
                        marker_code_found = part.split('=')[1].strip('[]')
                        break
            if marker_code_found == "N/A":
                marker_code_found = marker_code

            logging.info(f"  Processing record: ID='{record.id}'")
            logging.info(f"    Organism Found in Record: '{organism_name_found}'")
            logging.info(f"    Marker Found in Record (from header/fallback): '{marker_code_found}'")
            logging.info(f"    Checking for match with requested marker: '{marker_code}'")

            if not is_matching_marker(marker_code, marker_code_found):
                logging.warning(
                    f"  Skipped sequence '{record.id}' due to marker mismatch (final check): Found '{marker_code_found}' (canonical: {get_canonical_marker(marker_code_found)}), Expected '{marker_code}' (canonical: {get_canonical_marker(marker_code)})")
                continue

            canonical_marker = get_canonical_marker(marker_code_found) or marker_code_found

            sequences.append(
                SequenceData(
                    processid=processid,
                    organism_name=organism_name_found,
                    marker_code=canonical_marker,
                    sequence=str(record.seq),
                    query_name=organism
                )
            )
        logging.info(f"  Finished parsing all records from API response. Total sequences collected before limit: {len(sequences)}")

    except Exception as e:
        raise BOLDAPIError(f"Failed to parse FASTA response from BOLD Systems API. "
                           f"Error: {e}. "
                           f"Partial response: {response.text[:200]}...")

    if len(sequences) > limit:
        logging.info(f"Reducing {len(sequences)} found sequences to {limit} by random sampling.")
        sequences = random.sample(sequences, limit)
    else:
        logging.info(f"Number of sequences ({len(sequences)}) is within limit ({limit}). No sampling needed.")

    return sequences

if __name__ == "__main__":
    logging.info("--- BOLD Systems DNA Sequence Fetcher ---")

    organism_input = input("Enter organism names separated by commas (e.g. Homo sapiens, Mus musculus):\n> ")
    organisms_to_fetch = [org.strip() for org in organism_input.split(',') if org.strip()]

    common_barcode_type = input("Enter the barcode DNA type (e.g. COI-5P, 16S, ITS):\n> ").strip()
    if get_canonical_marker(common_barcode_type) is None:
        logging.warning(
            f"Entered marker '{common_barcode_type}' is not in the known alias list. Proceeding anyway, but be aware of potential mismatches.")

    try:
        sequences_per_organism_limit = int(input("Enter the number of sequences to fetch per organism:\n> ").strip())
        if sequences_per_organism_limit <= 0:
            raise ValueError("Limit must be positive.")
    except ValueError:
        logging.warning("Invalid number entered or limit must be positive. Using default value of 5.")
        sequences_per_organism_limit = 5

    all_fetched_sequences = []

    output_base_dir = os.path.join("data", "temp_fasta")

    canonical_marker_code = get_canonical_marker(common_barcode_type) or common_barcode_type
    combined_fasta_filename = f"combined_sequences_for_{canonical_marker_code.replace('-', '_')}.fasta"  # Replace hyphen for valid filename
    combined_fasta_path = os.path.join(output_base_dir, combined_fasta_filename)

    os.makedirs(output_base_dir, exist_ok=True)
    logging.info(f"Output directory ensured: {output_base_dir}")

    for organism in organisms_to_fetch:
        logging.info(f"Fetching {sequences_per_organism_limit} '{canonical_marker_code}' sequences for: {organism}")
        try:
            current_organism_sequences = fetch_sequences_by_organism_and_marker(
                organism=organism,
                marker_code=canonical_marker_code,
                limit=sequences_per_organism_limit
            )

            if current_organism_sequences:
                logging.info(f"  Successfully fetched {len(current_organism_sequences)} sequences for {organism}.")
                all_fetched_sequences.extend(current_organism_sequences)  # Add to the master list
            else:
                logging.warning(f"  No sequences found for {organism}.")
        except BOLDAPIError as e:
            logging.error(f"  Error fetching sequences for {organism}: {e}")
        except Exception as e:
            logging.exception(f"  An unexpected error occurred for {organism}.")

    if all_fetched_sequences:
        logging.info(f"--- Total sequences fetched: {len(all_fetched_sequences)} ---")
        logging.info(f"Attempting to create combined FASTA file at: {combined_fasta_path}")
        try:
            final_fasta_path = create_fasta_from_sequences(all_fetched_sequences, combined_fasta_path)
            logging.info(f"SUCCESS: All sequences collected into one FASTA file at: {final_fasta_path}")
        except SequenceProcessingError as e:
            logging.error(f"FAILURE: Error creating combined FASTA file: {e}")
        except Exception as e:
            logging.exception(f"FAILURE: An unexpected error occurred during FASTA file creation.")
    else:
        logging.warning("No sequences were fetched in total. No combined FASTA file will be created.")

    logging.info("--- Process Finished ---")
