"""
This program calculate the BBB score acccording to Gupa (2019) BBB Score

Usage: python BBB_Score_RS.py   (update the input file name in the code)
Input: input.sdf, a SD file, with a pKa_mb column.

Return: output.sdf, a SD file, with an added BBB_Score_RS column.
"""



import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as Descriptor
from rdkit.Chem import SDMolSupplier, SDWriter

def BBB_Descriptor(mol):
    """Calculate descriptors as per Gupa(2019) BBB Score"""
    # Calculate NWHBN
    nHBA = Descriptor.CalcNumHBA(mol)
    nHBD = Descriptor.CalcNumHBD(mol)
    HBN = nHBA + nHBD
    MW = round(Descriptor.CalcExactMolWt(mol, onlyHeavy=False), 2)
    MWHBN = round(HBN / (MW ** 0.5), 2)
    
    # Calculate number of the heavy atoms
    n = mol.GetNumAtoms()
    elements = [atom.GetSymbol() for atom in mol.GetAtoms()]
    nH = elements.count('H')
    HA = n - nH
    
    # Calculate the number of aromatic rings
    aroR = Descriptor.CalcNumAromaticRings(mol)
    
    # Calculate TPSA
    tpsa = Descriptor.CalcTPSA(mol)
    
    BBB_DESC = [MW, nHBA, nHBD, HBN, MWHBN, HA, aroR, tpsa]
    return BBB_DESC

def bbbscore(mol, pka):
    """Calculate BBB SCORE"""
    BBB_DESC = BBB_Descriptor(mol)
    HA = float(BBB_DESC[5])
    MWHBN = float(BBB_DESC[4])
    Aro_R = float(BBB_DESC[6])
    TPSA = float(BBB_DESC[7])
    
    # Calculate RS value for Aromatic Rings 
    if Aro_R == 0:
        RS_ARO_R = 0.336367
    elif Aro_R == 1:
        RS_ARO_R = 0.816016
    elif Aro_R == 2:
        RS_ARO_R = 1
    elif Aro_R == 3:
        RS_ARO_R = 0.691115
    elif Aro_R == 4:
        RS_ARO_R = 0.199399
    elif Aro_R > 4:
        RS_ARO_R = 0
    
    # Calculate RS value for HA
    if HA > 5 and HA <= 45:
        RS_HA = (0.0000443 * (HA ** 3) - 0.004556 * (HA ** 2) + 0.12775 * HA - 0.463) / 0.624231
    else:
        RS_HA = 0
    
    # Calculate RS value for MWHBN
    if MWHBN > 0.05 and MWHBN <= 0.45:
        RS_MWHBN = (26.733 * (MWHBN ** 3) - 31.495 * (MWHBN ** 2) + 9.5202 * MWHBN - 0.1358) / 0.72258
    else:
        RS_MWHBN = 0
    
    # Calculate RS value for TPSA
    if TPSA > 0 and TPSA <= 120:
        RS_TPSA = (-0.0067 * TPSA + 0.9598) / 0.9598
    else:
        RS_TPSA = 0
    
    # Calculate RS value for pKa
    pka = float(pka)
    if pka > 3 and pka <= 11:
        RS_PKA = (0.00045068 * (pka ** 4) - 0.016331 * (pka ** 3) + 0.18618 * (pka ** 2) - 0.71043 * pka + 0.8579) / 0.597488
    else:
        RS_PKA = 0
    
    BBB_score_RS = RS_ARO_R + RS_HA + 1.5 * RS_MWHBN + 2 * RS_TPSA + 0.5 * RS_PKA
    
    BBB_DESC.append(pka)
    BBB_DESC.append(round(BBB_score_RS, 2))
    
    return BBB_DESC

# Read input data from SD file
input_file = "input.sdf"
supplier = SDMolSupplier(input_file)

# Initialize output data
output_data = []

# Process each molecule in the input file
for mol in supplier:
    if mol is None:
        continue
    
    name = mol.GetProp("_Name")
    smiles = Chem.MolToSmiles(mol)
    
    if mol.HasProp("pKa_mb"):
        pka_mb = mol.GetProp("pKa_mb")
        bbb_score_rs = bbbscore(mol, pka_mb)
        
        # Append SMILES and ID to the output data along with the calculated scores
        output_data.append([name, smiles] + bbb_score_rs)

# Convert output data to DataFrame
output_df = pd.DataFrame(output_data, columns=['Name', 'SMILES', 'MW', 'nHBA', 'nHBD', 'HBN', 'MWHBN', 'HA', 'aroR', 'TPSA', 'pKa_mb', 'BBB_Score_RS'])

# Write output data to SD file
output_file = "output.sdf"
writer = SDWriter(output_file)

for index, row in output_df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    
    if mol is not None:
        mol.SetProp("_Name", row['Name'])
        mol.SetProp("MW", str(row['MW']))
        mol.SetProp("nHBA", str(row['nHBA']))
        mol.SetProp("nHBD", str(row['nHBD']))
        mol.SetProp("HBN", str(row['HBN']))
        mol.SetProp("MWHBN", str(row['MWHBN']))
        mol.SetProp("HA", str(row['HA']))
        mol.SetProp("aroR", str(row['aroR']))
        mol.SetProp("TPSA", str(row['TPSA']))
        mol.SetProp("pKa_mb", str(row['pKa_mb']))
        mol.SetProp("BBB_Score_RS", str(row['BBB_Score_RS']))
        
        writer.write(mol)

writer.close()
