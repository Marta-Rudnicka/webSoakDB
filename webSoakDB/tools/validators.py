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
from API.models import Compounds, Library, LibraryPlate, Project
from tools.uploads_downloads import standardize_smiles


#Library plate data uploads
def data_is_valid(file_name, error_log):
	'''Uses RDKit to validate SMILES strings, and includes RDKit
	error messages in the error log.'''
	
	print('entered data is valid')
	valid = True
	well_names = {}
	codes = {}
	
	if not is_csv(file_name, error_log):
		print('not a csv file')
		return False

	print('is a csv file')

	try:
		with open(file_name, newline='') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.read(1024))
			dialect.delimiter = ','
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
	inventory, and is it belongs to the selected library.'''
	
	valid = True
	
	if not is_csv(file_name, error_log):
		return False
	
	try:
		with open(file_name, newline='') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.read(1024))
			dialect.delimiter = ',' #because Sniffer gets confused with SMILES strings and requires manual correction
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
		code = row[0].strip()
		if not code_is_valid(code, line, codes, error_log):
			valid_row = False
	except (IndexError):
		return True 	#ignore empty lines
	
	try:
		well = row[1].strip()
		if not well_is_valid(well, line, well_names, error_log):
			valid_row = False
	except (IndexError):
		msg = "Line " + str(line) + ": CSV FORMATTING ERROR: Not enough fields. Line contains only 1 field, but 2 to 4 are required."
		update_error_log(msg, error_log)
		valid_row = False
		
	try:
		smiles = row[2].strip()
		if smiles != "" and smiles != " ":
			if not smiles_is_valid(smiles, line, error_log):
				valid_row = False
	except (IndexError):
		pass #this field is optional
	
	try:
		concentration = row[3].strip()
		if not concentration_is_valid(concentration, line, error_log):
			valid_row = False
	except (IndexError):
		pass	#this field is optional
	
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
		raw_smiles = row[0].strip()
	except (IndexError):
		return True # ignore empty lines

	if not smiles_is_valid(raw_smiles, line, error_log):
		return False

	smiles = standardize_smiles(raw_smiles)

	if not smiles_string_exists(smiles, line, error_log):
		valid_row = False
		
	try: 		#if a file has too many fields, it could be a wrong file that happen to have SMILES string in it
		extra_field = row[1]
		msg = "Line " + str(line) + ": CSV FORMATTING ERROR: Too many fields. Line contains more than 1 field."
		update_error_log(msg, error_log)
		valid_row = False
	except (IndexError):
		pass #everything is fine
	
	
	if not compound_exists(smiles, library_id, line, error_log):
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

def valid_well_name(string, line, error_log):

	if string == "":
		msg = "Line " + str(line) + ": WELL NAME ERROR. Missing well name."
		update_error_log(msg, error_log)
		return False

	match = re.fullmatch('([A-Z][A-F]?[1-9])|([A-Z][A-F]?(?!(00))([0-9]{2}))', string.upper())
	digit_part = '[0-9][0-9]?'
	
	if match and (int(re.search(digit_part, string)[0]) <= 48 ):
		return True
	
	msg = "Line " + str(line) + ": WELL NAME ERROR. Invalid well name: '" + string + "'"
	update_error_log(msg, error_log)
	return False

def unique_well_name(string, line, used_well_names, error_log):
	if string.upper() in used_well_names:
		msg = 'Line ' + str(line) + ": WELL NAME ERROR. Duplicate well name. Well name '" + string + "' was already used in line " + str(used_well_names[string.upper()])
		update_error_log(msg, error_log)
		used_well_names[string.upper()] = [used_well_names[string.upper()], line]
		return False
	
	used_well_names[string.upper()] = line
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
	'''If SMILES string is valid, returns empty string; otherwise returns RDKit error message'''
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
			c = float(string)
			if c > 0:
				return True
			else:
				msg = "Line {line}: CONCENTRATION ERROR. Negative concentration value ({string})".format(line=str(line), string=string)
				update_error_log(msg, error_log)
				return False
				
		except(ValueError):
			msg = "Line " + str(line) + ": CONCENTRATION ERROR. Concentration value '" + string + "' is not a number!"
			update_error_log(msg, error_log)
			return False
	return True

#FILE validation helpers

def is_csv(file_name, error_log):
	match = re.fullmatch('^.*\.csv$', file_name)

	if match:
		return True
	msg = "FILE ERROR: Wrong file type: '" + file_name + "'! Compound data should be uploaded as a CSV file."
	update_error_log(msg, error_log)
	return False

#COMPOUND finder


def compound_exists(smiles, library_id, line, error_log):
	if not smiles:
		return False

	matching_smiles = Compounds.objects.filter(smiles=smiles)

	if matching_smiles.count() == 0:
		msg = "Line " + str(line) + ": COMPOUND ERROR: '" + smiles + "': No such compound is registered in the inventory."
		update_error_log(msg, error_log)
		return False
	

	library_id = int(library_id)
	found = False

	for c in matching_smiles.all():
		if library_id in [sw.library_plate.library.id for sw in c.locations.all()]:
			found = True
		
	if not found:
		try:
			library = Library.objects.get(id=library_id)
		except django.core.exceptions.ObjectDoesNotExist:
			msg = "APPLICATION ERROR: Selected library does not exist"
			update_error_log(msg, error_log)
			return False
			
		msg = "Line " + str(line) + ": COMPOUND ERROR: no compound with the SMILES string: " + smiles + " belongs to " + library.name
		update_error_log(msg, error_log)
		return False
	return True



#EXPORT FORM FOR SOAKDB-COMPATIBLE LIST OF SOURCE COMPOUNDS

def export_form_is_valid(post_data):
	project = None
	subset_lib_ids = []
	for key in post_data:
		if key=='csrfmiddlewaretoken':
			pass
		elif key=='project':
			try:
				project = Project.objects.get(id=post_data.get(key))
				subset_lib_ids = [s.library.id for s in project.subsets.all()]
			except django.core.exceptions.ObjectDoesNotExist:
				print('No project found')
				return False
		else:
			if re.fullmatch('[0-9]+', key):
				value = post_data.get(str(key), False)
			elif re.fullmatch('([0-9]+)(\-[0-9]+)', key):
				value = post_data.get(str(key), False)
				old_key = re.fullmatch('([0-9]+)(\-[0-9]+)', key)
				key = old_key.group(1)
			else:
				print('key not matching any regex:', key)
				return False
				 
			if not int(key) in subset_lib_ids:
				print('key not found in subset libraries: ', key)
				return False
				
			try:
				plate = LibraryPlate.objects.get(pk=value)
			except django.core.exceptions.ObjectDoesNotExist:
				print('plate does not exist: ', value)
				return False
			
			if plate.library.id != int(key):
				print('plate library not matching subset library')
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
