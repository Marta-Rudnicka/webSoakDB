'''	Module to validate uploaded CSV files. Contains two validation scripts:
	* data_is_valid() - to validate compound data for a libary plate
	* selection_is_valid() - to validate cherrypicking list
	
	The scripts open the CSV file with file_name and validate the data 
	in the file. Whith each error, they add a new piece of HTML code with
	an error message to error_log. If the data is valid, they return True. 
	If errors have been found, or the file is not a CSV file, they return 
	False. They use csv.Sniffer().sniff() to detect the dialect of the 
	csv file.
	file_name: should be a string with the name of an existing file
	error_log: should be an empty list
		 
	Information about accepted data formats is at the end of the file.
'''

from rdkit import Chem
Chem.WrapLogs()
from io import StringIO
import sys
import re
import csv
import django.core.exceptions
from API.models import Compounds, Library


#Library plate data uploads
def data_is_valid(file_name, error_log):
	'''Uses RDKit to validate SMILES strings, and includes RDKit
	error messages in the error log.'''
	
	valid = True
	well_names = {}
	codes = {}
	
	if not is_csv(file_name, error_log):
		return False

	try:
		with open(file_name, newline='') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.read(1024))
			csvfile.seek(0)
			compound_reader = csv.reader(csvfile, dialect)
			line = 1
			for row in compound_reader:
				if not row_is_valid(row, line, error_log, well_names, codes):
					valid = False
				line += 1
		
	except FileNotFoundError:
		msg = "FILE ERROR: File '" + file_name + "' does not exist!"
		update_error_log(msg, error_log)
		return False
	
	return valid

#Cherrypicking list uploads
def selection_is_valid(file_name, error_log, library_id):
	'''Checks if file format is valid, if the compound exists in the
	inventory, and is it belongs to the selected library. Does not
	validate SMILES strings - the Compounds objects with invalid strings 
	should not exist in the DB anyway'''
	
	valid = True
	
	if not is_csv(file_name, error_log):
		return False
	
	try:
		with open(file_name, newline='') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.read(1024))
			csvfile.seek(0)
			compound_reader = csv.reader(csvfile, dialect)
			line = 1
			for row in compound_reader:
				if not compound_is_valid(row, line, error_log, library_id):
					valid = False
				line += 1
	except FileNotFoundError:
		msg = "FILE ERROR: File '" + file_name + "' does not exist!"
		update_error_log(msg, error_log)
		return False
	
	return valid

#EVERYTHING BELOW IS HELPERS#########################################

#for data_is_valid()
def row_is_valid(row, line, error_log, well_names, codes):
	valid_row = True
	
	try:
		code = row[0]
		if not code_is_valid(code, line, codes, error_log):
			valid_row = False
	except (IndexError):
		pass
	
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
		if smiles != "" and smiles != " ":
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


#for selection_is_valid()
def compound_is_valid(row, line, error_log, library_id):
	valid_row = True
	
	try:
		code = row[0].strip()
		if not code_exists(code, line, error_log):
			valid_row = False
	except (IndexError):
		pass #ignore empty lines
		
	try:
		smiles = row[1].strip()
		if not smiles_string_exists(smiles, line, error_log):
			valid_row = False
	except (IndexError):
		msg = "Line " + str(line) + ": CSV FORMATTING ERROR: Not enough fields. Line contains only 1 field, but 2 are required."
		update_error_log(msg, error_log)
		return False
	
	try:
		extra_field = row[2]
		msg = "Line " + str(line) + ": CSV FORMATTING ERROR: Too many fields. Line contains more than 2 fields."
		update_error_log(msg, error_log)
		valid_row = False
	except (IndexError):
		pass
	
	
	if not compound_exists(code, smiles, library_id, line, error_log):
		return False
	
	return valid_row

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


def smiles_string_exists(string, line, error_log):
	if string == "" or string == None:
		msg = 'Line ' + str(line) + ': SMILES STRING ERROR: Missing SMILES string.'
		update_error_log(msg, error_log)
		return False
	return True

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

#COMPOUND finder


def compound_exists(code, smiles, library_id, line, error_log):
	if not code or not smiles:
		return False

	try:
		compound = Compounds.objects.get(code=code, smiles=smiles)
	except django.core.exceptions.ObjectDoesNotExist:
		msg = "Line " + str(line) + ": COMPOUND ERROR: '" + code + ' : ' + smiles + "': No such compound is registered in the inventory."
		update_error_log(msg, error_log)
		return False
	except django.core.exceptions.MultipleObjectsReturned:
		msg = "Line " + str(line) + ": DATA INTEGRITY ERROR DETECTED: This selection cannot be uploaded because of database error. There are duplicate compounds in the database: ", row[0], ':', row[1], 'and this selection cannot be uploaded. Please contact the staff.'
	
	locations = compound.locations.all()
	library_id = int(library_id)
	found = False
	
	for location in locations:
		if location.library_plate.library.id == library_id:
			found = True
	
	if not found:
		try:
			library = Library.objects.get(id=library_id)
		except django.core.exceptions.ObjectDoesNotExist:
			msg = "Line " + str(line) + ": DATA ERROR: Selected library does not exist"
			update_error_log(msg, error_log)
			return False
			
		msg = "Line " + str(line) + ": COMPOUND ERROR: '" + code + ' : ' + smiles + "' does not belong to " + library.name
		update_error_log(msg, error_log)
		return False
	return True



'''
ACCEPTED DATA FORMATS:

Compound data for a library plate; accepted CSV data format:
	Field 1 (Code): any non-empty string
	Field 2 (Well name): one or two letters + one or two digits (except '0' and '00')
	Field 3 (SMILES string): valid SMILES string or empty string
	Field 4 (Concentration): integer, float, or a string that can be converted to a float, 
		empty string, or nothing(optional field)

Cherrypicking list: accepted CSV data format:
	Field 1 (Code): any non-empty string
	Field 2 (SMILES string): a non-empty string
'''
