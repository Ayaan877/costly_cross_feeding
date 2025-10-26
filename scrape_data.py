from Bio.KEGG import REST
import pandas as pd
import time
from tqdm import tqdm
import pickle

#-------------------------------------------------------------------------
# 1. Fetch reaction IDs from a KEGG pathway/module
#-------------------------------------------------------------------------
def get_reactions(pathway_id="map01100"):
    """
    Fetch all KEGG reaction IDs from a given pathway.
    Returns a list of reaction IDs.
    """
    response = REST.kegg_link("reaction", f"path:{pathway_id}").read()
    lines = response.strip().split("\n")
    reaction_ids = [line.split("\t")[1].split(":")[1] for line in lines]
    return reaction_ids

#-------------------------------------------------------------------------
# 2. Fetch the equation for a given reaction
#-------------------------------------------------------------------------
def get_equations(rxn_id):
    """
    Fetch the equation line for a KEGG reaction ID.
    Returns the equation string.
    """
    try:
        entry = REST.kegg_get(rxn_id).read()
        for line in entry.split("\n"):
            if line.startswith("EQUATION"):
                return line.replace("EQUATION", "").strip()
    except Exception as e:
        print(f"Failed to fetch {rxn_id}: {e}")
    return None

#-------------------------------------------------------------------------
# 3. Parse an equation into metabolites & coefficients
#-------------------------------------------------------------------------
def parse_equation(equation):
    """
    Parse a KEGG equation string into stoichiometry.
    Returns lists of (metabolites, coefficients).
    Reactants are negative, products are positive.
    Cancels out metabolites that appear on both sides.
    """
    stoich = {}

    if "<=>" in equation:
        lhs, rhs = equation.split(" <=> ")
    elif "=>" in equation:
        lhs, rhs = equation.split(" => ")
    else:
        return [], []

    # Reactants (negative)
    for term in lhs.split(" + "):
        parts = term.strip().split()
        if len(parts) == 0:
            continue
        coeff = int(parts[0]) if parts[0].isdigit() else 1
        met = parts[-1]
        if met.startswith(("C", "G")):
            stoich[met] = stoich.get(met, 0) - coeff

    # Products (positive)
    for term in rhs.split(" + "):
        parts = term.strip().split()
        if len(parts) == 0:
            continue
        coeff = int(parts[0]) if parts[0].isdigit() else 1
        met = parts[-1]
        if met.startswith(("C", "G")):
            stoich[met] = stoich.get(met, 0) + coeff

    # Drop any metabolites with net coefficient 0
    mets, coeffs = [], []
    for met, coeff in stoich.items():
        if coeff != 0:
            mets.append(met)
            coeffs.append(coeff)

    return mets, coeffs

#-------------------------------------------------------------------------
# 4. Build stoichiometric matrix from reactions
#-------------------------------------------------------------------------
def build_stoich_matrix(reaction_ids, sleep_time=0.1):
    """
    Build a stoichiometric matrix (metabolites x reactions).
    Each reversible reaction is split into forward and reverse columns.
    Returns a pandas DataFrame.
    """
    reaction_stoich = {}
    all_metabolites = set()
    
    for rxn_id in tqdm(reaction_ids, desc="Processing reactions"):
        eqn = get_equations(rxn_id)
        if not eqn:
            continue

        # Parse forward direction
        mets, coeffs = parse_equation(eqn)
        if mets:
            forward_id = f"{rxn_id}_f"
            reaction_stoich[forward_id] = dict(zip(mets, coeffs))
            all_metabolites.update(mets)

        # If reversible, also add reverse direction
        if "<=>" in eqn:
            reverse_id = f"{rxn_id}_r"
            reaction_stoich[reverse_id] = {m: -c for m, c in zip(mets, coeffs)}
            all_metabolites.update(mets)

        time.sleep(sleep_time)

    # Build DataFrame
    stoich_matrix = pd.DataFrame(
        0, index=sorted(all_metabolites), columns=reaction_stoich.keys()
    )
    for rxn, stoich in reaction_stoich.items():
        for met, coeff in stoich.items():
            stoich_matrix.loc[met, rxn] = coeff

    return stoich_matrix

#-------------------------------------------------------------------------
# 5. Scrape data and generate matrix
#-------------------------------------------------------------------------
def get_matrix(pathway_id):

    print("Fetching reactions...")
    rxns = get_reactions(pathway_id)
    print(f"Found {len(rxns)} reactions.")

    print("Building stoichiometric matrix...")
    stoich_matrix = build_stoich_matrix(rxns)

    print("Stoichiometric matrix shape:", stoich_matrix.shape)
    print(stoich_matrix.head())

    # Save to CSV
    stoich_matrix.to_csv(f"{pathway_id}_stoich_matrix.csv")
    print(f"Matrix saved to {pathway_id}_stoich_matrix.csv")

    return stoich_matrix

#-------------------------------------------------------------------------
# 6. Get compound and reaction strings from KEGG IDs
#-------------------------------------------------------------------------
def get_cpd_names(met_ids, sleep_time=0.1):
    """
    Builds a dictionary mapping KEGG compound IDs to their chemical names.
    """
    cpd_string_dict = {}
    for kegg_id in tqdm(met_ids, desc="Fetching compound names"):
        try:
            # Fetch entry from KEGG
            record = REST.kegg_get(kegg_id).read()
            # Parse the NAME field (first name is primary)
            for line in record.split('\n'):
                if line.startswith("NAME"):
                    name = line.split("NAME")[1].strip().split(';')[0]
                    cpd_string_dict[kegg_id] = name
                    break
            else:
                cpd_string_dict[kegg_id] = "Unknown"
        except Exception as e:
            cpd_string_dict[kegg_id] = "Unknown"
        time.sleep(sleep_time)  # Polite delay
    return cpd_string_dict

def get_rxn_names(rxn_ids, sleep_time=0.1):
    """
    Builds a dictionary mapping KEGG reaction IDs to their equation strings.
    """
    rxn_string_dict = {}
    base_eqn_cache = {}

    for rxn_id in tqdm(rxn_ids, desc="Fetching reaction equations"):
        # Remove _f or _r suffix for KEGG lookup
        base_id = rxn_id[:-2] if rxn_id.endswith(('_f', '_r')) else rxn_id

        # Only fetch from KEGG if not already cached
        if base_id not in base_eqn_cache:
            try:
                record = REST.kegg_get(base_id).read()
                equation = None
                for line in record.split('\n'):
                    if line.startswith("EQUATION"):
                        equation = line.split("EQUATION")[1].strip()
                        break
                if equation is None:
                    equation = "Unknown"
                base_eqn_cache[base_id] = equation
            except Exception as e:
                base_eqn_cache[base_id] = "Unknown"
            time.sleep(sleep_time)  # Polite delay

        rxn_string_dict[rxn_id] = base_eqn_cache[base_id]

    return rxn_string_dict

#-------------------------------------------------------------------------
# 7. Run script
#-------------------------------------------------------------------------
if __name__ == "__main__":
    print("Choose an option:")
    print("1. Generate stoichiometric matrix from KEGG pathway")
    print("2. Generate compound and reaction dictionaries from existing matrix CSV")
    print("3. Do nothing / Exit")
    choice = input("Enter 1, 2 or 3: ").strip()

    if choice == "1":
        pathway_id = input("Enter KEGG pathway ID (e.g., map01100): ").strip()
        stoich_matrix = get_matrix(pathway_id)
    elif choice == "2":
        pathway_id = input("Enter KEGG pathway ID (e.g., map01100): ").strip()
        csv_file = f"{pathway_id}_stoich_matrix.csv"
        df = pd.read_csv(csv_file, index_col=0)
        met_ids = df.index.tolist()
        rxn_ids = df.columns.tolist()

        print("Generating compound name dictionary...")
        cpd_string_dict = get_cpd_names(met_ids)
        with open("cpd_string_dict.pkl", "wb") as f:
            pickle.dump(cpd_string_dict, f)
        print("Saved cpd_string_dict.pkl")

        print("Generating reaction equation dictionary...")
        rxn_string_dict = get_rxn_names(rxn_ids)
        with open("rxn_string_dict.pkl", "wb") as f:
            pickle.dump(rxn_string_dict, f)
        print("Saved rxn_string_dict.pkl")
    elif choice == "3":
        print("No action taken. Exiting.")
    else:
        print("Invalid choice.")
