from django.test import TestCase
from API.models import Library, LibraryPlate, Compounds, SourceWell
import unittest
from tools.compounds import find_compound_in_plate
from tools.set_up import (
    set_up_source_wells,
	source_wells_data,
	plates_data,
	libraries_data,
	compounds_data,
)

class TestFindCompoundInPlate(unittest.TestCase):
    def setUp(self):
        set_up_source_wells(source_wells_data, plates_data, libraries_data, compounds_data)
    
    def test1(self):
        lib1 = Library.objects.get(name="lib1")
        plate1 = lib1.plates.all()[0]
        c1 = Compounds.objects.filter(code="code2")[0]
        c2 = Compounds.objects.filter(code="code9")[0]
        sw1 = SourceWell.objects.get(library_plate=plate1, compound = c1)

        self.assertEqual(find_compound_in_plate(c1, plate1), sw1)
        self.assertEqual(find_compound_in_plate(c2, plate1), None)

    def test_alternatives(self):
        #lib1 = Library.objects.get(name="lib1")
        plate1 = LibraryPlate.objects.get(barcode="xyz")
        plate2 = LibraryPlate.objects.get(barcode="xyz_v2")
        c1 = Compounds.objects.filter(code="code1")[0]
        c2 = Compounds.objects.filter(code="code2")[0]
        c3 = Compounds.objects.filter(code="code8")[0]
        
        self.assertEqual(find_compound_in_plate(c1, plate1).compound.code, "code1")
        self.assertEqual(find_compound_in_plate(c1, plate2).compound.code, "alt_code1")
        self.assertEqual(find_compound_in_plate(c2, plate1).compound.code, "code2")
        self.assertEqual(find_compound_in_plate(c3, plate2), None) #no alternative compound