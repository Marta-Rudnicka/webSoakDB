from API.models import Library, LibraryPlate, Compounds, SourceWell, Proposals, LibrarySubset
import django.db
import csv
import django.core.exceptions

#import data from a csv file - should be used after validation with data_is_valid
def get_dict_from_echo_file(file_name, plate):
	with open(file_name, newline='') as csvfile:
		dialect = csv.Sniffer().sniff(csvfile.read(1024))
		csvfile.seek(0)
		compound_reader = csv.reader(csvfile, dialect)
		print(compound_reader[0])
		#for row in compound_reader:
			
