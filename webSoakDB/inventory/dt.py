'''scripts used for dispense testing'''
from webSoakDB_backend.validators import is_csv
from API.models import SWStatuschange, SourceWell
import csv
import string
import re
from django.core.exceptions import ObjectDoesNotExist

def create_well_dict():
	
	#'human readable' format: [A - H][01-12][a-d], e.g. B09c, F12a
	
	#generate lists of strings to concatenate into well names
	part_1 = list(string.ascii_uppercase)[0:8]
	part_2 = []
	for x in range(1,13):
		if x<10:
			part_2.append('0' + str(x))
		else:
			part_2.append(str(x))

	part_3 = [['a', 'c'], ['b', 'd']]

	#generate an ordered list of all well names in the human readable system
	well_name_list_1 = []
	
	for number in part_2:
		for charlist in part_3:
			for cap in part_1:
				for letter in charlist:
					well_name_list_1.append(cap + number + letter)

	
	#Echo  format: [A-P][1-24], e.g D2, H19
	
	#generate an ordered list of all well names
	well_name_list_echo = []

	part_1_echo = list(string.ascii_uppercase)[0:16]
	for number in range(1, 25):
		for letter in part_1_echo:
			well_name_list_echo.append(letter + str(number))


	#create a dictionary matching the names
	well_dictionary = {}

	i = 0
	for item in well_name_list_1:
		well_dictionary[well_name_list_echo[i]] = item
		i = i + 1

	return well_dictionary

def ensure_leading_zero(string):
	'''in human-readable format, e.g. turn A4d to A04d'''
	pattern = '([A-H])([1-9])([a,c,d])'
	if re.fullmatch(pattern, string):
		name = re.fullmatch(pattern, string)
		return name.group(1) + '0' + name.group(2) + name.group(3)
	else:
		return string

#validate well names
def valid_wells(dw, sw):
	
	pattern_echo = '[A-P][1-9][0-9]?'
	pattern_hr = '[A-H][0-9][0-9]?[a,c,d]'
	pattern_sw = '[A-Z][A-F]?[0-9]{1,2}'
	
	digit_part = '[0-9][0-9]?'
	
	error_string = ""
	if not (re.fullmatch(pattern_echo, dw) or re.fullmatch(pattern_hr, dw)):
		error_string = 'Invalid destination well: ' + dw + "; "
	
	if not re.fullmatch(pattern_sw, sw):
		error_string = error_string + 'Invalid source well: ' + sw + "; "
	
	if (re.fullmatch(pattern_echo, dw) and int(re.search(digit_part, dw)[0]) > 24 ):
		error_string = 'Invalid destination well: ' + dw + "; "
	
	if (re.fullmatch(pattern_hr, dw) and int(re.search(digit_part, dw)[0]) > 12 ):
		error_string = 'Invalid destination well: ' + dw + "; "
	
	if (re.fullmatch(pattern_sw, sw) and int(re.search(digit_part, sw)[0]) > 48 ):
		error_string = error_string + 'Invalid source well: ' + sw + "; "
		
	if error_string == "":
		return True
	else:
		return error_string
	
#map destination wells to source wells based on input file
def get_well_dictionary(file_name):
	
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		reader = csv.reader(csvfile, dialect)
		
		#find columns with source and destination wells, return error messages if not found
		header = [s.strip().lower() for s in next(reader, None)]
		try:
			sw = header.index('source well')
		except(ValueError):
			return "Column for source well not found. Make sure your file uses the header 'Source well' for your source well column."
		try:
			dw = header.index('destination well')
		except(ValueError):
			return "Column for destination well not found. Make sure your file uses the header 'Destination well' for your destination well column"
					
		#determine naming system for destination well based on the first destination well name
		name = next(reader, None)[dw]
		
		pattern_echo = '[A-P][1-9][0-9]?'
		pattern_hr = '[A-H][0-9][0-9]?[a,c,d]'
		
		if re.fullmatch(pattern_echo, name):
			formatting = "echo"
			well_name_converter = create_well_dict()
		elif re.fullmatch(pattern_hr, name):
			formatting = "human-readable"
		else:
			return "Unrecognised formatting of the name of the destination well (based on the first well)."	
		
		well_dict = {}
		error_log = []
		
		#go back to the second line of the file (first line of data)
		csvfile.seek(0)
		next(reader)
		
		#map destination wells to source wells, or create an error log
		for row in reader:
			if row[dw]=="" and row[sw] == "":
				continue
			if valid_wells(row[dw], row[sw]) == True:
				try:
					if formatting == "echo" and re.fullmatch(pattern_echo, row[dw]):
							well_dict[well_name_converter[row[dw]]] = row[sw]
					elif formatting == "human-readable" and re.fullmatch(pattern_hr, row[dw]):
						well_dict[ensure_leading_zero(row[dw])] = row[sw]
					else: 
						error_log.append('Unrecognised well: ' + row[dw] + '(does not fit the pattern from the first well)')
				
				except (KeyError):
					error_log.append('Unrecognised well: ' + row[dw] + '(does not fit the pattern from the first well)')
			
			else:
				error_log.append(valid_wells(row[dw], row[sw]))
		
		if error_log == []:
			return well_dict
		else:
			return error_log

def manage_sw_status_change(source_well, date, activation):
	try:
		change = SWStatuschange.objects.get(source_well = source_well, date=date)
		#reverse to the state before the change
		if change.activation != activation:
			source_well.active = activation
			change.delete()
			if activation:
				source_well.deactivation_date = None
			else:
				last_change = source_well.status_changes.order_by("date").reverse()[0]
				source_well.deactivation_date = last_change.date
			
			source_well.save()
		return
	except(ObjectDoesNotExist):
		SWStatuschange.objects.create(source_well=source_well, date=date, activation=activation)
		if activation:
			source_well.deactivation_date = None
		else:
			source_well.deactivation_date = date
		source_well.save()
