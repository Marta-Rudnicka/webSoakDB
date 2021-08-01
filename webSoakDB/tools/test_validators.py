from django.test import TestCase
import unittest
import validators as v

from API.models import Library, SourceWell

from tools.set_up import (
    set_up_source_wells,
	source_wells_data,
	plates_data,
	libraries_data,
	compounds_data,
)


class TestRowIsValid(unittest.TestCase):
    
    def test_valid_input(self):
        error_log = []
        well_names = {}
        codes = {}

        valid_inputs = [
            ["code_string1", "A4", "C(NC1CC1)C2CCOC2", "20"],
            ["code_string2", "A5", "C(NC1CC1)C2CCOC2" ],
            ["code_string3", "A6"],
            ["code_string4", "A7", ""],
            ["code_string5", "A8", "", ""],
            ["code_string6", "A9", "", "98"],
            []
        ]

        i = 1
        for row in valid_inputs:
            self.assertTrue(v.row_is_valid(row, i, error_log, well_names, codes))
            i += 1

        self.assertEqual(error_log, [])
        final_codes = {
            "code_string1" : 1, 
            "code_string2" : 2, 
            "code_string3" : 3, 
            "code_string4" : 4, 
            "code_string5" : 5,
            "code_string6" : 6
        }
        self.assertEqual(codes, final_codes)
    
    def test_missing_values(self):
        self.maxDiff = None
        error_log = []
        well_names = {}
        codes = {}

        inputs = [
            ["", "A4", "C(NC1CC1)C2CCOC2", "20"],
            ["", "A5", "", "20"],
            ["", "F45", "", ""],
            ["", "G02", ""],
            ["", "G06"],
            ["code_string1", "", "C(NC1CC1)C2CCOC2", "80" ],
            ["code_string2", ""],
            ["code_string3", "", ""],
            ["code_string4", "", "CC", ""],
            ["code_string5", "", "CC"],
            ["code_string6", "", "", "90"],
            ["", ""]
        ]

        i = 1
        for row in inputs:
            self.assertFalse(v.row_is_valid(row, i, error_log, well_names, codes))
            i += 1

        errs = [
            "<p>Line 1: COMPOUND CODE ERROR: Missing compound code.</p>",
            "<p>Line 2: COMPOUND CODE ERROR: Missing compound code.</p>",
            "<p>Line 3: COMPOUND CODE ERROR: Missing compound code.</p>",
            "<p>Line 4: COMPOUND CODE ERROR: Missing compound code.</p>",
            "<p>Line 5: COMPOUND CODE ERROR: Missing compound code.</p>",
            "<p>Line 6: WELL NAME ERROR. Missing well name.</p>",
            "<p>Line 7: WELL NAME ERROR. Missing well name.</p>",
            "<p>Line 8: WELL NAME ERROR. Missing well name.</p>",
            "<p>Line 9: WELL NAME ERROR. Missing well name.</p>",
            "<p>Line 10: WELL NAME ERROR. Missing well name.</p>",
            "<p>Line 11: WELL NAME ERROR. Missing well name.</p>",
            "<p>Line 12: COMPOUND CODE ERROR: Missing compound code.</p>",
            "<p>Line 12: WELL NAME ERROR. Missing well name.</p>",
        ]

        #print(error_log)
        self.assertEqual(error_log, errs)
        
    def test_wrong_formatting(self):
        self.maxDiff = None
        error_log = []
        well_names = {}
        codes = {}

        inputs = [
            ["code_string1", "A4", "C(NC1CC1)C2CCOC2", "20", ""],
            ["code_string2", "A5", "CC", "", ""],
            ["code_string3", "F45", "", "", ""],
            ["code_string4", "A8", "", "20", ""],
            ["code_string5"]
        ]

        i = 1
        for row in inputs:
            self.assertFalse(v.row_is_valid(row, i, error_log, well_names, codes))
            i += 1

        errs = [
            "<p>Line 1: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
            "<p>Line 2: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
            "<p>Line 3: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
            "<p>Line 4: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
            "<p>Line 5: CSV FORMATTING ERROR: Not enough fields. Line contains only 1 field, but 2 to 4 are required.</p>"
        ]

        self.assertEqual(error_log, errs)

    def test_invalid_values(self):

            self.maxDiff = None
            error_log = []
            well_names = {}
            codes = {}

            inputs = [
                ["", "A4", "C(NC1CC1)C2CCOC2", "20"],
                ["code_string2", "A057", "CC", "80"],
                ["code_string3", "F42", "string", ""],
                ["code_string4", "A8", "CC", "string"],
                
            ]

            i = 1
            for row in inputs:
                self.assertFalse(v.row_is_valid(row, i, error_log, well_names, codes))
                i += 1

            errs = [
                "<p>Line 1: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
                "<p>Line 2: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
                "<p>Line 3: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
                "<p>Line 4: CSV FORMATTING ERROR: Too many fields. Line contains more than 4 fields.</p>",
            ]

            #self.assertEqual(error_log, errs)

    def test_duplicate_values(self):
        self.maxDiff = None
        error_log = []
        well_names = {}
        codes = {}

        valid_inputs = [
            ["code_string1", "A4"],
            ["code_string2", "A05"]                
        ]

        invalid_inputs = [
            ["code_string1", "F42" ],
            ["code_string3", "A4"],
            ["code_string1", "A05"],
                
        ]

        i = 1
        for row in valid_inputs:
            self.assertTrue(v.row_is_valid(row, i, error_log, well_names, codes))
            i += 1
        for row in invalid_inputs:
            self.assertFalse(v.row_is_valid(row, i, error_log, well_names, codes))
            i += 1
        errs = [
            "<p>Line 3: COMPOUND CODE ERROR: Duplicate. Code 'code_string1' was already used in line 1</p>",
            "<p>Line 4: WELL NAME ERROR. Duplicate well name. Well name 'A4' was already used in line 1</p>",
            "<p>Line 5: COMPOUND CODE ERROR: Duplicate. Code 'code_string1' was already used in line [1, 3]</p>",
            "<p>Line 5: WELL NAME ERROR. Duplicate well name. Well name 'A05' was already used in line 2</p>",
            
        ]

        self.assertEqual(error_log, errs)

class TestCompoundIsValid(unittest.TestCase):
    def setUp(self):
        set_up_source_wells(source_wells_data, plates_data, libraries_data, compounds_data)
    
    def test_valid_inputs(self):
        error_log = []
        library_id1 = Library.objects.get(name="lib1").id
        library_id2 = Library.objects.get(name="lib3").id
        
        rows1 = [
            ["CC"],
            ["CCO"],
            ["CCI"]
        ]

        rows2 = [
            ["CCC"],
            ["CCCF"],
            ["CCOF"]
        ]

        i = 1
        for row in rows1:
            self.assertTrue(v.compound_is_valid(row, i, error_log, library_id1))
            i += 1
        
        for row in rows2:
            self.assertTrue(v.compound_is_valid(row, i, error_log, library_id2))
            i += 1
        self.assertEqual(error_log, [])
    
    def test_invalid_inputs(self):
        print('test_invalid_inputs')
        error_log = []
        library_id1 = Library.objects.get(name="lib1").id

        rows = [
            ["CCCF"],
            ["string"],
            [""],
            ["CC", "string"]
        ]

        i = 1
        for row in rows:
            self.assertFalse(v.compound_is_valid(row, i, error_log, library_id1))
            i += 1

        errs = [
            "<p>Line 1: COMPOUND ERROR: no compound with the SMILES string: CCCF belongs to lib1</p>",
            "rdkit error - do not test",
            '<p>Line 3: SMILES STRING ERROR: Missing SMILES string.</p>',
            '<p>Line 4: CSV FORMATTING ERROR: Too many fields. Line contains more than 1 field.</p>'
        ]

        self.assertEqual(error_log[0], errs[0])
        self.assertEqual(error_log[2], errs[2])
        self.assertEqual(error_log[3], errs[3])

class TestWellValidation(TestCase):
    def test_valid_well_name(self):
        error_log = []

        valid_wells = ['A1', 'Z1', 'A48', 'Z48', 'AF1', 'AF48', 'B01', 'AD01']
        for str in valid_wells:
            self.assertTrue(v.valid_well_name(str, 1, error_log))
        
        invalid_wells = ['', 'A0', 'B00', 'F123', 'AG01', 'AB123', 'F49']
        i = 1

        for str in invalid_wells:
            self.assertFalse(v.valid_well_name(str, i, error_log))
            i += 1
        
        errs = [
           "<p>Line 1: WELL NAME ERROR. Missing well name.</p>",
           "<p>Line 2: WELL NAME ERROR. Invalid well name: 'A0'</p>",
           "<p>Line 3: WELL NAME ERROR. Invalid well name: 'B00'</p>",
           "<p>Line 4: WELL NAME ERROR. Invalid well name: 'F123'</p>",
           "<p>Line 5: WELL NAME ERROR. Invalid well name: 'AG01'</p>",
           "<p>Line 6: WELL NAME ERROR. Invalid well name: 'AB123'</p>",
           "<p>Line 7: WELL NAME ERROR. Invalid well name: 'F49'</p>"
       ]

        self.assertEqual(error_log, errs)
    
    def test_unique_well_name(self):
        self.maxDiff = None
        error_log = []
        used_well_names = {}
        wells = ['A1', 'Z1', 'A48', 'Z48', 'AF1', 'AF48', 'B01', 'AD01']
        
        i = 1
        for str in wells:
            self.assertTrue(v.unique_well_name(str, i, used_well_names, error_log))
            i += 1

        for str in wells:
            self.assertFalse(v.unique_well_name(str, i, used_well_names, error_log))
            i += 1
        
        errs = [
            "<p>Line 9: WELL NAME ERROR. Duplicate well name. Well name 'A1' was already used in line 1</p>",
            "<p>Line 10: WELL NAME ERROR. Duplicate well name. Well name 'Z1' was already used in line 2</p>",
            "<p>Line 11: WELL NAME ERROR. Duplicate well name. Well name 'A48' was already used in line 3</p>",
            "<p>Line 12: WELL NAME ERROR. Duplicate well name. Well name 'Z48' was already used in line 4</p>",
            "<p>Line 13: WELL NAME ERROR. Duplicate well name. Well name 'AF1' was already used in line 5</p>",
            "<p>Line 14: WELL NAME ERROR. Duplicate well name. Well name 'AF48' was already used in line 6</p>",
            "<p>Line 15: WELL NAME ERROR. Duplicate well name. Well name 'B01' was already used in line 7</p>",
            "<p>Line 16: WELL NAME ERROR. Duplicate well name. Well name 'AD01' was already used in line 8</p>",
        ]

        self.assertEqual(error_log, errs)

    def test_well_is_valid(self):
        error_log = []
        used_well_names = {}
        valid_wells = ['A1', 'Z1', 'A48', 'Z48', 'AF1', 'AF48', 'B01', 'AD01']

        i = 1
        for str in valid_wells:
            self.assertTrue(v.well_is_valid(str, i, used_well_names, error_log))
            i += 1
        
        invalid_wells = ['', 'A0', 'B00', 'F123', 'AG01', 'AB123', 'F49', 'A1', 'A0']
        for str in invalid_wells:
            self.assertFalse(v.well_is_valid(str, i, used_well_names, error_log))
            i += 1
        
        errs = [
           "<p>Line 9: WELL NAME ERROR. Missing well name.</p>",
           "<p>Line 10: WELL NAME ERROR. Invalid well name: 'A0'</p>",
           "<p>Line 11: WELL NAME ERROR. Invalid well name: 'B00'</p>",
           "<p>Line 12: WELL NAME ERROR. Invalid well name: 'F123'</p>",
           "<p>Line 13: WELL NAME ERROR. Invalid well name: 'AG01'</p>",
           "<p>Line 14: WELL NAME ERROR. Invalid well name: 'AB123'</p>",
           "<p>Line 15: WELL NAME ERROR. Invalid well name: 'F49'</p>",
           "<p>Line 16: WELL NAME ERROR. Duplicate well name. Well name 'A1' was already used in line 1</p>",
           "<p>Line 17: WELL NAME ERROR. Invalid well name: 'A0'</p>",
        ]

        #print(error_log)
        self.assertEqual(error_log, errs)
        
class TestSmilesValidation(TestCase):

    def test_smiles_string_exists(self):
        error_log = []
        self.assertTrue(v.smiles_string_exists("string", 1, error_log))
        self.assertFalse(v.smiles_string_exists("", 2, error_log))
        self.assertFalse(v.smiles_string_exists(None, 3, error_log))

        errs = [
            '<p>Line 2: SMILES STRING ERROR: Missing SMILES string.</p>',
            '<p>Line 3: SMILES STRING ERROR: Missing SMILES string.</p>',
        ]

        self.assertEqual(error_log, errs)


    def test_parse_smiles(self):
        valid1 = "C(NC1CC1)C2CCOC2"
        valid2 = " C(NC1CC1)C2CCOC2 "
        non_smiles = "badstring"
        invalid = "CN(C)(C)C"
        self.assertEqual(v.parse_smiles(valid1), '')
        self.assertEqual(v.parse_smiles(valid2), '')
        self.assertNotEqual(v.parse_smiles(non_smiles), '')
        self.assertNotEqual(v.parse_smiles(invalid), '')

    def test_smiles_is_valid(self):
        self.maxDiff = None
        error_log = []
        valid1 = "C(NC1CC1)C2CCOC2"
        valid2 = " C(NC1CC1)C2CCOC2 "
        non_smiles = "badstring"
        invalid = "CN(C)(C)C"

        self.assertTrue(v.smiles_is_valid(valid1, 1, error_log))
        self.assertTrue(v.smiles_is_valid(valid2, 2, error_log))
        self.assertFalse(v.smiles_is_valid(non_smiles, 3, error_log))
        self.assertFalse(v.smiles_is_valid(invalid, 4, error_log))

        self.assertAlmostEqual(len(error_log), 2)
        #Error messages contain timestamps, so it's impossible to predict the correct contents

class TestConcentrationIsValid(TestCase):
    def test(self):
        error_log = []
        self.assertTrue(v.concentration_is_valid("100", 1, error_log))
        self.assertTrue(v.concentration_is_valid("3.15", 2, error_log))
        self.assertTrue(v.concentration_is_valid("0100", 3, error_log))
        self.assertTrue(v.concentration_is_valid(456, 4, error_log))
        self.assertTrue(v.concentration_is_valid("", 5, error_log))
        self.assertFalse(v.concentration_is_valid("str", 6, error_log))
        self.assertFalse(v.concentration_is_valid("-40", 7, error_log))
        self.assertFalse(v.concentration_is_valid(-40, 8, error_log))

        errs = [
            "<p>Line 6: CONCENTRATION ERROR. Concentration value 'str' is not a number!</p>",
            "<p>Line 7: CONCENTRATION ERROR. Negative concentration value (-40)</p>",
            "<p>Line 8: CONCENTRATION ERROR. Negative concentration value (-40)</p>"
        ]

        self.assertEqual(error_log, errs)

class TestIsCSV(TestCase):
    def test(self):
        error_log = []
        self.assertTrue(v.is_csv("x.csv", error_log))
        self.assertTrue(v.is_csv("x-n- 0.csv", error_log))
        self.assertTrue(v.is_csv("x.xls.csv", error_log))
        self.assertFalse(v.is_csv("x_- i .pdf", error_log))
        self.assertFalse(v.is_csv("x_- i .csv.xls", error_log))

        errs = [
            "<p>FILE ERROR: Wrong file type: 'x_- i .pdf'! Compound data should be uploaded as a CSV file.</p>",
            "<p>FILE ERROR: Wrong file type: 'x_- i .csv.xls'! Compound data should be uploaded as a CSV file.</p>",
        ]

        self.assertEqual(error_log, errs)

class TestCompoundExists(TestCase):
    def setUp(self):
        set_up_source_wells(source_wells_data, plates_data, libraries_data, compounds_data)
    
    def test_valid(self):
        l = Library.objects.get(name="lib1")
        l_smiles = ["CC", "CCO", "CCI", "CCF"]
        error_log = []
        i=1
        for str in l_smiles:
            self.assertTrue(v.compound_exists(str, l.id, i, error_log))
            i += 1

    def test_wrong_lib(self):
        self.maxDiff = None
        l = Library.objects.get(name="lib2")
        smiles = ["CC", "CCO", "CCI", "CCF"]
        error_log = []
        i=1
        for str in smiles:
            self.assertFalse(v.compound_exists(str, l.id, i, error_log))
            i += 1

        errs = [
            "<p>Line 1: COMPOUND ERROR: no compound with the SMILES string: CC belongs to lib2</p>",
            "<p>Line 2: COMPOUND ERROR: no compound with the SMILES string: CCO belongs to lib2</p>",
            "<p>Line 3: COMPOUND ERROR: no compound with the SMILES string: CCI belongs to lib2</p>",
            "<p>Line 4: COMPOUND ERROR: no compound with the SMILES string: CCF belongs to lib2</p>"
        ]
    
        self.assertEqual(error_log, errs)
    

    def test_smiles_not_in_db(self):
        self.maxDiff = None
        l = Library.objects.get(name="lib1")
        smiles = ["str1", "str2"]
        error_log = []
        i=1
        for str in smiles:
            self.assertFalse(v.compound_exists(str, l.id, i, error_log))
            i += 1
        
        errs = [
            "<p>Line 1: COMPOUND ERROR: 'str1': No such compound is registered in the inventory.</p>",
            "<p>Line 2: COMPOUND ERROR: 'str2': No such compound is registered in the inventory.</p>"
        ]
    
        self.assertEqual(error_log, errs)
    
    def test_app_error(self):
        error_log =[]
        self.assertFalse(v.compound_exists("CCO", 99, 1, error_log))
        err = [
            "<p>APPLICATION ERROR: Selected library does not exist</p>"
        ]
        self.assertEqual(error_log, err)

class TestCodeValidation(TestCase):
    def test_code_exists(self):
        error_log = []
        
        self.assertFalse(v.code_exists("", 1, error_log))
        self.assertTrue(v.code_exists("str", 2, error_log))

        errs = ["<p>Line 1: COMPOUND CODE ERROR: Missing compound code.</p>"]

        self.assertEqual(error_log, errs)
    
    def test_unique_code(self):
        self.maxDiff = None
        error_log = []
        used_codes = {}
        codes = ['code1', 'code2']
        
        i = 1
        for str in codes:
            self.assertTrue(v.unique_code(str, i, used_codes, error_log))
            i += 1

        for str in codes:
            self.assertFalse(v.unique_code(str, i, used_codes, error_log))
            i += 1
        
        errs = [
            "<p>Line 3: COMPOUND CODE ERROR: Duplicate. Code 'code1' was already used in line 1</p>",
            "<p>Line 4: COMPOUND CODE ERROR: Duplicate. Code 'code2' was already used in line 2</p>",
        ]

        self.assertEqual(error_log, errs)
    
    def test_code_is_valid(self):
        self.maxDiff = None
        error_log = []
        used_codes = {}
        valid_codes = ['code1', 'code2']
        invalid_codes = ['code1', 'code2', ""]
        
        i = 1
        for str in valid_codes:
            self.assertTrue(v.code_is_valid(str, i, used_codes, error_log))
            i += 1

        for str in invalid_codes:
            self.assertFalse(v.code_is_valid(str, i, used_codes, error_log))
            i += 1
        
        errs = [
            "<p>Line 3: COMPOUND CODE ERROR: Duplicate. Code 'code1' was already used in line 1</p>",
            "<p>Line 4: COMPOUND CODE ERROR: Duplicate. Code 'code2' was already used in line 2</p>",
            "<p>Line 5: COMPOUND CODE ERROR: Missing compound code.</p>"
        ]

        self.assertEqual(error_log, errs)
