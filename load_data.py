import pandas as pd
import numpy as np
import pickle

'''
Build dictionaries for internal reference IDs and KEGG IDs.
Returns:
    stoich_matrix: 2D stoichiometric matrix array

    met_map: KEGG metabolite ID to internal index
    inv_met_map: internal index to KEGG metabolite ID
    rxn_map: KEGG reaction ID to internal index
    inv_rxn_map: internal index to KEGG reaction ID
    cpd_string_dict: KEGG metabolite IDs to chemical compound names
    display_lookup: KEGG reaction IDs to reaction strings

    mets: list of internal indices for all metabolites
    rxns: list of internal indices for all reactions
    Currency: list of internal indices for currency metabolites
    Energy: list of internal indices for nutrient metabolites
    Core: list of internal indices for core metabolites

    rho: reactant stoichiometric matrix (negative entries of stoich_matrix, others set to 0)
    pi: product stoichiometric matrix (positive entries of stoich_matrix, others set to 0)
    rxnMat: binary reactant matrix (1 if metabolite is a reactant in reaction, else 0)
    prodMat: binary product matrix (1 if metabolite is a product in reaction, else 0)
    sumRxnVec: vector of counts of reactants per reaction
    sumProdVec: vector of counts of products per reaction
'''

# Load stoichiometric matrix
stoich_matrix_df = pd.read_csv("map01100_stoich_matrix.csv", index_col=0)
stoich_matrix = stoich_matrix_df.values

metabolites = stoich_matrix_df.index.tolist()
reactions = stoich_matrix_df.columns.tolist()

# Metabolite dictionaries
met_map = kegg_to_id = {met_id: idx for idx, met_id in enumerate(metabolites)}
inv_met_map = id_to_kegg = {idx: met_id for met_id, idx in met_map.items()}
with open(f"met_map.pkl", "wb") as f:
    pickle.dump(met_map, f)
with open(f"inv_met_map.pkl", "wb") as f:
    pickle.dump(inv_met_map, f)

# Reactions dictionaries
rxn_map = rxn_kegg_to_id = {rxn_id: idx for idx, rxn_id in enumerate(reactions)}
inv_rxn_map = rxn_id_to_kegg = {idx: rxn_id for rxn_id, idx in rxn_map.items()}
with open(f"rxn_map.pkl", "wb") as f:
    pickle.dump(rxn_map, f)
with open(f"inv_rxn_map.pkl", "wb") as f:
    pickle.dump(inv_rxn_map, f)

# Currency and nutrient metabolites
indices_of_currency_mets = pd.read_csv("kegg_currency.txt", header=None)[0].tolist()
indices_of_energy_mets = pd.read_csv("kegg_nutrients.txt", header=None)[0].tolist()
indices_of_core_mets = pd.read_csv("kegg_core.txt", header=None)[0].tolist()

Currency = [met_map[i] for i in indices_of_currency_mets if i in met_map]
Energy = [met_map[i] for i in indices_of_energy_mets if i in met_map]
Core = [met_map[i] for i in indices_of_core_mets if i in met_map]

# Compound names dictionary
mets = list(met_map.values())
rxns = list(rxn_map.values())
with open("cpd_string_dict.pkl", "rb") as f:
    cpd_string_dict = pickle.load(f)
with open("rxn_string_dict.pkl", "rb") as f:
    display_lookup = pickle.load(f)

# Defining reactant, product and reactant vectors/matrices for the scope expansion algorithm.
rho = stoich_matrix.clip(max = 0.0)
pi = stoich_matrix.clip(min = 0.0)
rxnMat = (rho != 0) * 1
prodMat = (pi != 0) * 1
sumRxnVec = np.sum(rxnMat, axis = 1)
sumProdVec = np.sum(prodMat, axis = 1)

# #-------------------------------------------------------------------------
# # Check for non-reversible reactions
# #-------------------------------------------------------------------------

# from collections import defaultdict

# # Group reactions by their base KEGG ID (remove suffix)
# rxn_groups = defaultdict(list)
# for rxn_id in reactions:  # "reactions" = stoich_matrix.columns.tolist()
#     base_id = rxn_id[:-2]  # strip "_f" or "_r"
#     rxn_groups[base_id].append(rxn_id)

# # Keep only those that have a single direction (irreversible)
# non_reversible_rxns = [rxns[0] for rxns in rxn_groups.values() if len(rxns) == 1]

# # Get their internal indices
# non_reversible_indices = [rxn_map[rxn_id] for rxn_id in non_reversible_rxns]

# print("Number of non-reversible reactions:", len(non_reversible_rxns))
# print("Example KEGG IDs:", non_reversible_rxns)
