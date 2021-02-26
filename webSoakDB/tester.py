#from API.models import Library, Compounds
from webSoakDB_backend.validators import selection_is_valid
from temp_code_storage.dupl import duplicate_codes

'''
def run(file_name, library_id):
	#for debugging
	error_log = []
	value = selection_is_valid(file_name, error_log, library_id)
	
	print('return value of is_valid: ', value)
	print('eror log: ')
	for msg in error_log:
		print(msg)
'''

def run(file_name):
	duplicate_codes(file_name)
