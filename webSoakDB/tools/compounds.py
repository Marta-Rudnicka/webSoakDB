import django.core.exceptions
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

def standardize_smiles(raw_string):
	'''Canonicalize SMILES and remove salts - use only on validated SMILES strings;
	For parsing CSV files'''

	#remove whitespace and convert to canonical form
	smiles = Chem.CanonSmiles(raw_string.strip())

	#convert to RDKit mol object 
	mol = Chem.MolFromSmiles(smiles)

	#remove salt
	remover = SaltRemover()
	desalted_mol = remover.StripMol(mol, dontRemoveEverything=True)

	#get SMILES string of the desalted mol and return it
	return Chem.MolToSmiles(desalted_mol)

def find_compound_in_plate(compound, library_plate):
    '''returns SourceWell in <library_plate> with <compound> in it ; if not found, 
    looks for a different Compounds object with the same SMILES string; if still 
    not found, returns None'''

    try:
        sw = library_plate.compounds.filter(compound = compound) #filter not to crash in case of duplicated data
        return sw[0]
    except IndexError:
        sw = find_alternative_compound(compound, library_plate)
        return sw

def find_alternative_compound(compound, plate):
    '''returns SourceWell with with a compound that has the same SMILES string
    as <compound>; if not found, returns None'''
    sws = plate.compounds.filter(compound__smiles=compound.smiles)
    if sws.count() > 0:
        return sws[0]
    else:
        return None

