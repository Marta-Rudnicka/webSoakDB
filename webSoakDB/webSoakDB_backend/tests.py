from django.test import TestCase
from inventory.set_up import (set_up_source_wells, source_wells_data, plates_data, libraries_data, compounds_data)
from tools.uploads_downloads import import_full_libraries
from API.models import Project, SourceWell, LibraryPlate

from tools.validators import code_is_valid, well_is_valid, smiles_is_valid, concentration_is_valid, row_is_valid

class ValidatorHelpersTests(TestCase):
	
	def test1(self):
		self.assertTrue(code_is_valid('randomstring', 2, {}, [] ))
	
	def test2(self):
		self.assertTrue(code_is_valid('123456879', 2, {}, [] ))

	def test3(self):
		self.assertTrue(code_is_valid("str", 2, {'abc' : 1, 'def' : 2}, [] ))
	
	def test4(self):
		self.assertFalse(code_is_valid('abc', 2, {'abc' : 1, 'def' : 2}, [] ))
	
	def test5(self):
		self.assertFalse(code_is_valid('', 2, {'abc' : 1, 'def' : 2}, [] ))
		
	def test6(self):
		self.assertTrue(well_is_valid("A01", 3, {}, []))

	def test7(self):
		self.assertTrue(well_is_valid("a1", 3, {}, []))

	def test8(self):
		self.assertTrue(well_is_valid("AF48", 3, {}, []))

	def test9(self):
		self.assertTrue(well_is_valid("b09", 3, {}, []))

	def test10(self):
		self.assertFalse(well_is_valid("abc13", 3, {}, []))

	def test11(self):
		self.assertFalse(well_is_valid("c18", 3, {"C18": 2}, []))

	def test12(self):
		self.assertFalse(well_is_valid("G49", 3, {}, []))

	def test13(self):
		self.assertFalse(well_is_valid("AG01", 3, {}, []))
	
	def test14(self):
		self.assertFalse(smiles_is_valid('randomstring', 2, []))

	def test15(self):
		self.assertFalse(smiles_is_valid('CCC=(O)O', 2, []))

	def test16(self):
		self.assertFalse(smiles_is_valid('68', 2, []))
	
	def test17(self):
		self.assertTrue(smiles_is_valid('', 2, []))

	def test18(self):
		self.assertTrue(smiles_is_valid('CC(=O)NCCC1=CNc2c1cc(OC)cc2', 2, []))

	def test19(self):
		self.assertTrue(concentration_is_valid('12', 2, []))
	
	def test20(self):
		self.assertTrue(concentration_is_valid('38.75', 2, []))
	
	def test21(self):
		self.assertTrue(concentration_is_valid('', 2, []))
	
	def test22(self):
		self.assertTrue(concentration_is_valid('015', 2, []))

	def test23(self):
		self.assertFalse(concentration_is_valid('-8', 2, []))
	
	def test24(self):
		self.assertFalse(concentration_is_valid('string', 2, []))

	def test25(self):
		self.assertTrue(concentration_is_valid(2014, 2, []))
	
	def test26(self):
		self.assertTrue(row_is_valid(['string', 'A3', 'CC', '100'], 3, [], {}, {}))

	def test27(self):
		self.assertFalse(row_is_valid(['', 'A3', 'CC', '100'], 3, [], {}, {}))
	
	def test28(self):
		self.assertFalse(row_is_valid(['string', 'AG3', 'foo', '100'], 3, [], {}, {}))
	
	def test29(self):
		self.assertFalse(row_is_valid(['string', 'A3', 'CC', 'bar'], 3, [], {}, {}))

class HelpersTests(TestCase):
	def setUp():
		#set_up()
		set_up_source_wells(source_wells_data, plates_data, libraries_data, compounds_data)
		
	def import_full_libraries_test():
		proposal = Project.objects.get(proposal="proposal1")
		p1 = LibraryPlate.objects.get(barcode="xyz")
		p2 = LibraryPlate.objects.get(barcode="abc1")
		
		sw1 = SourceWell.objects.get(library_plate=p1, well="A01")
		sw2 = SourceWell.objects.get(library_plate=p1, well="B01")
		sw3 = SourceWell.objects.get(library_plate=p1, well="C01")
		sw4 = SourceWell.objects.get(library_plate=p1, well="D01")
		sw5 = SourceWell.objects.get(library_plate=p4, well="A01")
		sw6 = SourceWell.objects.get(library_plate=p4, well="B01")
		sw7 = SourceWell.objects.get(library_plate=p4, well="C01")
		
		output = [sw1, sw2, sw3, sw4, sw5, sw6, sw7]
		
		self.assertEqual(import_full_libraries(proposal), output)

class ViewsTests(TestCase):
	pass

