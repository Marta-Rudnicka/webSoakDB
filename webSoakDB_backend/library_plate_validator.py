'''	Module to validate library plate data in a CSV file
	data_is_valid() performs the validation; the rest of the functions are helpers.
	Accepted CSV format:
	Field 1 (Code): any non-empty string
	Field 2 (Well name): one or two letters + one or two digits (except '0' and '00')
	Field 3 (SMILES string): valid SMILES string or empty string
	Field 4 (Concentration): integer, float, or a string that can be converted to a float, 
		empty string, or nothing(optional field)
	
'''

from rdkit import Chem
from io import StringIO
import sys
import re
import csv
Chem.WrapLogs()


def data_is_valid(file_name, error_log):
	'''Opens CSV file with file_name and validates the data in the file.
	Whith each error it finds, it adds a new string to error_log 
	(the strings are HTML code). If the data is valid, returns True.
	If errors have been found, or the file is not a CSV file, returns
	False. Uses RDKit to validate SMILES strings, and includes RDKit
	error messages in the error log. Input: file_name should be a string, and
	error_log should be a list (preferable an empty one). Uses Python
	csv.Sniffer().sniff() function to detect the format (dialect) of the csv file.'''
	
	valid = True
	well_names = {}
	codes = {}
	
	if not is_csv(file_name, error_log):
		valid = False
		return valid

	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		line = 1
		for row in compound_reader:
			if not row_is_valid(row, line, error_log, well_names, codes):
				valid = False
			line += 1
	return valid

def run(file_name):
	'''for debugging'''
	error_log = []
	value = is_valid(file_name, error_log)
	
	print('return value of is_valid: ', value)
	print('eror log: ')
	for msg in error_log:
		print(msg)


def update_error_log(error_message, log):
	line = "<p>" + error_message + "</p>"
	log.append(line)

#CODE validation helpers

def code_exists(string, line, error_log):
	if string == "":
		msg = 'Line ' + str(line) + ': COMPOUND CODE ERROR: Missing compound code.'
		update_error_log(msg, error_log)
		return False
	return True

def unique_code(string, line, used_codes, error_log):
	if string in used_codes:
		msg = 'Line ' + str(line) + ": COMPOUND CODE ERROR: Duplicate. Code '" + string + "' was already used in line " + str(used_codes[string])
		update_error_log(msg, error_log)
		used_codes[string] = [used_codes[string], line]
		return False
	used_codes[string] = line
	return True

def code_is_valid(string, line, used_codes, error_log):
	return code_exists(string, line, error_log) and unique_code(string, line, used_codes, error_log)

#WELL validation helpers

def valid_well_name(string, line, error_log): #figure out what are actual possible well names
	match = re.fullmatch('([A-Z]{1,2}[1-9])|([A-Z]{1,2}(?!(00))([0-9]{2}))', string.upper())
	if match:
		return True
	
	msg = "Line " + str(line) + ": WELL NAME ERROR. Invalid well name: '" + string + "'"
	update_error_log(msg, error_log)
	return False

def unique_well_name(string, line, used_well_names, error_log):
	if string in used_well_names:
		msg = 'Line ' + str(line) + ": WELL NAME ERROR. Duplicate well name. Well name '" + string + "' was already used in line " + str(used_well_names[string])
		update_error_log(msg, error_log)
		used_well_names[string] = [used_well_names[string], line]
		return False
	used_well_names[string] = line
	return True

def well_is_valid(string, line, used_well_names, error_log):
	return valid_well_name(string, line, error_log) and unique_well_name(string, line, used_well_names, error_log)


#SMILES validation helpers

def parse_smiles(string):
	sio = sys.stderr = StringIO()
	mol = Chem.MolFromSmiles(string)
	output = sio.getvalue()
	return output

def smiles_is_valid(string, line, error_log):
	rdkit_stderr = parse_smiles(string)
	if rdkit_stderr == "" and Chem.MolFromSmiles(string) != None:
		return True
	
	msg = "Line " + str(line) + ": SMILES STRING ERROR. Invalid SMILES string: '" + string + '\'<br/><code class="rdkit-err">' + rdkit_stderr + "</code>"
	update_error_log(msg, error_log)
	return False


#CONCENTRATION validation helpers

def concentration_is_valid(string, line, error_log):
	if string != "":
		try:
			float(string)
			return True
		except(ValueError):
			msg = "Line " + str(line) + ": CONCENTRATION ERROR. Concentration value '" + string + "' is not a number!"
			update_error_log(msg, error_log)
			return False
	return True

#FILE validation helpers

def is_csv(file_name, error_log):
	match = re.fullmatch('(.*)+\.csv$', file_name)

	if match:
		return True
	msg = "FILE ERROR: Wrong file type: '" + file_name + "'! Compound data should be uploaded as a CSV file."
	update_error_log(msg, error_log)
	return False

#validator helper

def row_is_valid(row, line, error_log, well_names, codes):
	valid_row = True
	
	try:
		code = row[0]
		if not code_is_valid(code, line, codes, error_log):
			valid_row = False
	except (IndexError):
		pass #ignore empty lines
	
	try:
		well = row[1]
		if not well_is_valid(well, line, well_names, error_log):
			valid_row = False
	except (IndexError):
		msg = "Line: " + str(line) + "CSV FORMATTING ERROR: Not enough fields. Line contains only 1 field, but 3 or 4 are required."
		update_error_log(msg, error_log)
		valid_row = False
		
	try:
		smiles = row[2]
		if smiles != "":
			if not smiles_is_valid(smiles, line, error_log):
				valid_row = False
	except (IndexError):
		msg = "Line " + str(line) + ": CSV FORMATTING ERROR: Not enough fields. Line contains only 2 fields, but 3 or 4 are required. If not providing the SMILES string, leave an empty field"
		update_error_log(msg, error_log)
		valid_row = False
	
	try:
		concentration = row[3]
		if not concentration_is_valid(concentration, line, error_log):
			valid_row = False
	except (IndexError):
		pass
	
	try:
		extra_field = row[4]
		msg = "Line " + str(line) + ": CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields."
		update_error_log(msg, error_log)
		valid_row = False
	except (IndexError):
		pass
	
	return valid_row
