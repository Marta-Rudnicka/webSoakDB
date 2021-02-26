'''	Module to validate cherrypicking lists submitted as a CSV file
	selection_is_valid() performs the validation; the rest of the functions are helpers.
	Accepted CSV format:
	Field 1 (Code): any non-empty string
	Field 2 (SMILES string): a non-empty string
	
'''

#from rdkit import Chem
#from io import StringIO
#import sys
import re
import csv
import django.core.exceptions
from API.models import Compounds, Library
#Chem.WrapLogs()


def selection_is_valid(file_name, error_log, library_id):
	'''Opens CSV file with file_name and validates the data in the file.
	Whith each error it finds, it adds a new string to error_log 
	(the strings are HTML code). If the data is valid, returns True.
	If errors have been found, or the file is not a CSV file, returns
	False. Input: file_name should be a string, error_log should be a list
	(preferable an empty one). Uses Python csv.Sniffer().sniff() function 
	to detect the format (dialect) of the csv file.'''
	
	valid = True
	
	if not is_csv(file_name, error_log):
		valid = False
		return valid

	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		line = 1
		for row in compound_reader:
			if not compound_is_valid(row, line, error_log, library_id):
				valid = False
			line += 1
	return valid

def run(file_name, library_id):
	'''for debugging'''
	error_log = []
	value = selection_is_valid(file_name, error_log, library_id)
	
	print('return value of is_valid: ', value)
	print('eror log: ')
	for msg in error_log:
		print(msg)


def update_error_log(error_message, log):
	line = "<p>" + error_message + "</p>"
	log.append(line)

#CODE validation helpers

def code_exists(string, line, error_log):
	if string == "" or string == None:
		msg = 'Line ' + str(line) + ': COMPOUND CODE ERROR: Missing compound code.'
		update_error_log(msg, error_log)
		return False
	return True

#######################################################			NEW
#SMILES validation helpers

def smiles_string_exists(string, line, error_log):
	if string == "" or string == None:
		msg = 'Line ' + str(line) + ': SMILES STRING ERROR: Missing SMILES string.'
		update_error_log(msg, error_log)
		return False
	return True
##########################################################		END NEW
#FILE validation helpers

def is_csv(file_name, error_log):
	match = re.fullmatch('(.*)+\.csv$', file_name)

	if match:
		return True
	msg = "FILE ERROR: Wrong file type: '" + file_name + "'! Compound data should be uploaded as a CSV file."
	update_error_log(msg, error_log)
	return False

#validator helper

def compound_is_valid(row, line, error_log, library_id):
	valid_row = True
	
	try:
		code = row[0]
		if not code_exists(code, line, error_log):
			valid_row = False
	except (IndexError):
		pass #ignore empty lines
		
	try:
		smiles = row[1]
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

def compound_exists(code, smiles, library_id, line, error_log):
	if not code or not smiles:
		return False

	try:
		compound = Compounds.objects.get(code=code, smiles=smiles)
	except django.core.exceptions.ObjectDoesNotExist:
		msg = "Line " + str(line) + ": COMPOUND ERROR: '" + code + ' : ' + smiles + "': No such compound is registered in the inventory."
		update_error_log(msg, error_log)
		return False
	
	locations = compound.locations.all()
	found = False
	for location in locations:
		if location.library_plate.library.id == library_id:
			found = True
	
	if not found:
		library = Library.objects.get(id=library_id)
		msg = "Line " + str(line) + ": COMPOUND ERROR: '" + code + ' : ' + smiles + "' does not belong to " + library.name
		update_error_log(msg, error_log)
		return False
	return True
