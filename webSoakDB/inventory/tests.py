from django.test import TestCase, Client
from django.contrib.auth.models import User
from django.urls import reverse
import string
import unittest
import datetime

from .dt import ensure_leading_zero, valid_wells
from API.models import Library, LibraryPlate, Preset, Compounds, SourceWell, SWStatuschange, PlateOpening, Proposals
from .inv_helpers import get_plate_size, fake_compounds_copy, get_change_dates, set_status, get_usage_stats, current_library_selection, fake_preset_copy

def set_up():
	l1 = Library.objects.create(name='lib1', public=True, for_industry=True)
	l2 = Library.objects.create(name='lib2', public=False, for_industry=True)
	l3 = Library.objects.create(name='lib3', public=True, for_industry=False)
	
	p1 = LibraryPlate.objects.create(barcode="xyz", library=l1, current=True)
	p2 = LibraryPlate.objects.create(barcode="xyz2", library=l1, current=False)
	p3 = LibraryPlate.objects.create(barcode="xyz3", library=l1, current=False)
	p4 = LibraryPlate.objects.create(barcode="abc1", library=l2, current=True)
	p5 = LibraryPlate.objects.create(barcode="largest_small", library=l2, current=False)
	p6 = LibraryPlate.objects.create(barcode="smallest_large", library=l2, current=False)
	p7 = LibraryPlate.objects.create(barcode="large_row", library=l2, current=False)
	p8 = LibraryPlate.objects.create(barcode="large_column", library=l2, current=False)
	p9 = LibraryPlate.objects.create(barcode="empty", library=l2, current=False)
	
		
	c1 = Compounds.objects.create(code="code1", smiles="CC")
	c2 = Compounds.objects.create(code="code2", smiles="CCO")
	c3 = Compounds.objects.create(code="code3", smiles="CCI")
	c4 = Compounds.objects.create(code="code4", smiles="CCF")
	c5 = Compounds.objects.create(code="code5", smiles="")
	c6 = Compounds.objects.create(code="code6", smiles="")
	c7 = Compounds.objects.create(code="code7", smiles="")
	c8 = Compounds.objects.create(code="code8", smiles="CCFO")
		
	sw1 = SourceWell.objects.create(compound=c1, library_plate=p1, well="A01", concentration=30, active=True)
	sw2 = SourceWell.objects.create(compound=c2, library_plate=p1, well="B01", concentration=30, active=True)
	sw3 = SourceWell.objects.create(compound=c3, library_plate=p1, well="C01", concentration=30, active=True)
	sw4 = SourceWell.objects.create(compound=c4, library_plate=p1, well="D01", concentration=30, active=True)
	sw5 = SourceWell.objects.create(compound=c1, library_plate=p2, well="A01", concentration=30, active=False, deactivation_date=datetime.date(2021, 2, 10))
	sw6 = SourceWell.objects.create(compound=c2, library_plate=p2, well="B01", concentration=30, active=True)
	sw7 = SourceWell.objects.create(compound=c3, library_plate=p2, well="C01", concentration=30, active=False, deactivation_date=datetime.date(2020, 10, 28))
	sw8 = SourceWell.objects.create(compound=c4, library_plate=p2, well="D01", concentration=30, active=True)
	sw9 = SourceWell.objects.create(compound=c1, library_plate=p3, well="A01", concentration=30, active=True)
	sw10 = SourceWell.objects.create(compound=c2, library_plate=p3, well="B01", concentration=30, active=True)
	sw11 = SourceWell.objects.create(compound=c3, library_plate=p3, well="C01", concentration=30, active=False, deactivation_date=datetime.date(2020, 10, 28))
	sw12 = SourceWell.objects.create(compound=c4, library_plate=p3, well="D01", concentration=30, active=True)
	sw13 = SourceWell.objects.create(compound=c5, library_plate=p4, well="A01", concentration=30, active=True)
	sw14 = SourceWell.objects.create(compound=c6, library_plate=p4, well="B01", concentration=30, active=True)
	sw15 = SourceWell.objects.create(compound=c7, library_plate=p4, well="C01", concentration=30, active=False, deactivation_date=datetime.date(2020, 10, 28))
	
	sw14 = SourceWell.objects.create(compound=c8, library_plate=p5, well="P24", concentration=30, active=True)
	sw15 = SourceWell.objects.create(compound=c8, library_plate=p6, well="R25", concentration=30, active=True)
	sw16 = SourceWell.objects.create(compound=c8, library_plate=p7, well="R02", concentration=30, active=True)
	sw17 = SourceWell.objects.create(compound=c8, library_plate=p8, well="F25", concentration=30, active=True)
	
	sch1 = SWStatuschange.objects.create(source_well=sw7, activation=False, date=datetime.date(2020, 10, 28))
	sch2 = SWStatuschange.objects.create(source_well=sw11, activation=False, date=datetime.date(2020, 10, 28))
	sch3 = SWStatuschange.objects.create(source_well=sw5, activation=False, date=datetime.date(2020, 10, 28))
	sch4 = SWStatuschange.objects.create(source_well=sw5, activation=True, date=datetime.date(2020, 12, 2))
	sch5 = SWStatuschange.objects.create(source_well=sw5, activation=False, date=datetime.date(2021, 2, 10))
		
	po1 = PlateOpening.objects.create(plate=p2, date="2020-10-28", reason="reason")
	po2 = PlateOpening.objects.create(plate=p3, date="2020-10-28", reason="reason")
	
	pr = Preset.objects.create(name="preset1", description="desc")
		
	prop = Proposals.objects.create(proposal="proposal1")
	prop.libraries.add(l1)
	prop.libraries.add(l2)
	prop.save()

class DTHelperTestsNoDB(unittest.TestCase):
	'''helper functions for dispense testing that do not use db objects)'''
	def test1(self):
		self.assertTrue(valid_wells("A01a", "A1"))
	
	def test2(self):
		self.assertTrue(valid_wells("H12d", "AF48"))
	
	def test3(self):
		self.assertTrue(valid_wells("A1", "A01"))

	def test4(self):
		self.assertTrue(valid_wells("P24", "F34"))
	
	def test5(self):
		self.assertEqual(valid_wells("A04b", "B8"), "Invalid destination well: A04b; ")
	
	def test6(self):
		self.assertEqual(valid_wells("I06c", "AG38"), "Invalid destination well: I06c; Invalid source well: AG38; ")
	
	def test7(self):
		self.assertEqual(valid_wells("A13b", "B49"), "Invalid destination well: A13b; Invalid source well: B49; ")
	
	def test8(self):
		self.assertEqual(valid_wells("P25", "Z88"), "Invalid destination well: P25; Invalid source well: Z88; ")
	
	def test9(self):
		self.assertEqual(valid_wells("PS5", "zjk7"), "Invalid destination well: PS5; Invalid source well: zjk7; ")
	
	def test10(self):
		self.assertEqual(ensure_leading_zero('A1d'), 'A01d')

	def test11(self):
		self.assertEqual(ensure_leading_zero('A03d'), 'A03d')
		
	def test12(self):
		self.assertEqual(ensure_leading_zero('H11d'), 'H11d')

	def test13(self):
		self.assertEqual(ensure_leading_zero('foo'), 'foo')
	
class HelperTestsDB(TestCase):
	'''inventory helper functions using DB objects'''
	def setUp(self):
		set_up()

	def fake_compounds_copy_test(self):
		l = LibraryPlate.objects.get(barcode="xyz2")
		queryset = l.compounds.all()
		copy = fake_compounds_copy(queryset)
		
		changes0 = [{"date" : datetime.date(2020, 10, 28), "activation" : False}, {"date" : datetime.date(2020, 12, 2), "activation" : True}, {"date" : datetime.date(2021, 2, 10), "activation" : False}]
		changes1 = []
		changes2 = [{"date" : datetime.date(2020, 10, 28), "activation" : False}]
		changes3 = []
		
		self.assertEqual(copy[0].well, "A01" )
		self.assertEqual(copy[1].well, "B01" )
		self.assertEqual(copy[2].well, "C01" )
		self.assertEqual(copy[3].well, "D01" )
				
		self.assertEqual(copy[0].code, "code1")
		self.assertEqual(copy[1].code, "code2")
		self.assertEqual(copy[2].code, "code3")
		self.assertEqual(copy[3].code, "code4")
		
		self.assertFalse(copy[0].active)
		self.assertTrue(copy[1].active)
		self.assertFalse(copy[2].active)
		self.assertTrue(copy[3].active)
		
		self.assertEqual(copy[0].deactivation_date, datetime.date(2021, 2, 10))
		self.assertEqual(copy[1].deactivation_date, None)
		self.assertEqual(copy[2].deactivation_date, datetime.date(2020, 10, 28))
		self.assertEqual(copy[3].deactivation_date, None)
		
		self.assertEqual(copy[0].changes, changes0)
		self.assertEqual(copy[1].changes, changes1)
		self.assertEqual(copy[2].changes, changes2)
		self.assertEqual(copy[3].changes, changes3)
	
	def fake_preset_copy_test(self):
		pass
			
	def get_plate_size_test(self):
		'''get_plate_size'''
		plate1 = LibraryPlate.objects.get(barcode='xyz')
		plate2 = LibraryPlate.objects.get(barcode='largest_small')
		plate3 = LibraryPlate.objects.get(barcode='smallest_large')
		plate4 = LibraryPlate.objects.get(barcode='large_row')
		plate5 = LibraryPlate.objects.get(barcode='large_column')
		
		queryset1 = SourceWell.objects.filter(library_plate=plate1)
		queryset2 = plate2.compounds.all()
		queryset3 = plate3.compounds.all()
		queryset4 = plate4.compounds.all()
		queryset5 = plate5.compounds.all()
		
		self.maxDiff = None
		
		small_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']
		small_columns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']
		small_plate = {'rows' : small_rows, 'columns' : small_columns}
		large_rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA', 'AB', 'AC','AD', 'AE', 'AF']
		large_columns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48']
		large_plate = {'rows' : large_rows, 'columns' : large_columns}
		self.assertEqual(get_plate_size(queryset1), small_plate)
		self.assertEqual(get_plate_size(queryset2), small_plate)
		self.assertEqual(get_plate_size(queryset3), large_plate)
		self.assertEqual(get_plate_size(queryset4), large_plate)
		self.assertEqual(get_plate_size(queryset5), large_plate)
		
	def get_change_dates_test(self):
		
		date201028 = datetime.date(2020, 10, 28)
		date201202 = datetime.date(2020, 12, 2)
		date210210 = datetime.date(2021, 2, 10)
		
		plate1 = LibraryPlate.objects.get(barcode="xyz2")
		compounds1 = plate1.compounds.all()
		input1 = fake_compounds_copy(compounds1)
		dates1 = [date210210, date201202, date201028]
	
		plate2 = LibraryPlate.objects.get(barcode="largest_small")
		compounds2 = plate2.compounds.all()
		input2 = fake_compounds_copy(compounds2)
		dates2 = []
		
		plate3 = LibraryPlate.objects.get(barcode="xyz3")
		compounds3 = plate3.compounds.all()
		input3 = fake_compounds_copy(compounds3)
		dates3 = [date201028]

		self.assertEqual(get_change_dates(input1), dates1)
		self.assertEqual(get_change_dates(input2), dates2)
		self.assertEqual(get_change_dates(input3), dates3)
	
	def set_status_test(self):
		'''set_status'''
		plate = LibraryPlate.objects.get(barcode="xyz2")
		sw = SourceWell.objects.get(library_plate=plate, well="A01")
		sw_copy = fake_compounds_copy([sw])
		
		date1 = datetime.date(2020, 10, 27)
		date2 = datetime.date(2020, 10, 28)#
		date3 = datetime.date(2020, 12, 1)
		date4 = datetime.date(2020, 12, 2)#
		date5 = datetime.date(2020, 12, 3)
		date6 = datetime.date(2021, 2, 9)
		date7 = datetime.date(2021, 2, 10)#
		date8 = datetime.date(2021, 2, 11)
		
		output1 = set_status(sw_copy, date1)
		output2 = set_status(sw_copy, date2)
		output3 = set_status(sw_copy, date3)
		output4 = set_status(sw_copy, date4)
		output5 = set_status(sw_copy, date5)
		output6 = set_status(sw_copy, date6)
		output7 = set_status(sw_copy, date7)
		output8 = set_status(sw_copy, date8)
		
		self.assertTrue(output1.active)
		self.assertFalse(output2.active)
		self.assertFalse(output3.active)
		self.assertTrue(output4.active)
		self.assertTrue(output6.active)
		self.assertFalse(output7.active)
		self.assertFalse(output8.active)

	def get_usage_stats_test(self):
		'''get_usage_stats'''
		plate1 = LibraryPlate.objects.get(barcode="xyz2")
		compounds1 = fake_compounds_copy(plate1.compounds.all())
		stats1 = {"count" : 4, "active" : 2, "inactive" : 2, "availability" : 50 }
		
		plate2 = LibraryPlate.objects.get(barcode="xyz3")
		compounds2 = fake_compounds_copy(plate2.compounds.all())
		stats2 = {"count" : 4, "active" : 3, "inactive" : 1, "availability" : 75 }
		
		plate3 = LibraryPlate.objects.get(barcode="empty")
		compounds3 = fake_compounds_copy(plate3.compounds.all())
		stats3 = {"count" : 0, "active" : 0, "inactive" : 0, "availability" : 0}
		
		self.assertEqual(get_usage_stats(compounds1), stats1)
		self.assertEqual(get_usage_stats(compounds2), stats2)
		self.assertEqual(get_usage_stats(compounds3), stats3)
	
	def current_library_selection_test(self):
		'''current_library_selection'''
		output_true = [("", "Select library..."), (1, "lib1"), (2, "lib2"), (3, "lib3")]
		output_false = [(1, "lib1"), (2, "lib2"), (3, "lib3")]
		self.assertEqual(current_library_selection(True), output_true)
		self.assertEqual(current_library_selection(False), output_false)
			
class GetViewsTests(TestCase):
	
	def setUp(self):
		set_up()
		self.user = User.objects.create_user('staff', email='staff@example.com', password='password', is_staff=True)
			
	def test1(self):
		'''inventory index'''
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('index'))
		self.assertEqual(response.status_code, 200)
		
	def test2(self):
		'''libraries'''
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('libraries'))
		self.assertEqual(response.status_code, 200)
		self.assertEqual(response.context["libraries"].count(), 2)
	
	def test3(self):
		'''lib. deletion error'''
		c = Client()
		response = c.get("/inventory/library-deletion-error/")
		self.assertEqual(response.status_code, 200)
	
	def test4(self):
		'''plates'''
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('plates'))
		self.assertEqual(response.status_code, 200)
		self.assertEqual(response.context["libraries"].count(), 2)
	
	def test5(self):
		'''presets'''
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('presets'))
		self.assertEqual(response.status_code, 200)
		#preset count
	
	def test7(self):
		'''plate lookup - full plate'''
		l1 = Library.objects.get(name="lib1")
		l = LibraryPlate.objects.get(barcode="xyz", library=l1)
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('update-plate', args=[l.id]))
		self.assertEqual(response.status_code, 200)
		self.assertEqual(response.context["compounds"].count(), 4)
		self.assertEqual(response.context["active_count"], 4)
		self.assertEqual(response.context["inactive_count"], 0)
		self.assertEqual(response.context["availability"], 100)

	def test8(self):
		'''plate lookup - partially used'''
		l1 = Library.objects.get(name="lib1")
		l = LibraryPlate.objects.get(barcode="xyz2", library=l1)
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('update-plate', args=[l.id]))
		self.assertEqual(response.status_code, 200)
		self.assertEqual(response.context["compounds"].count(), 4)
		self.assertEqual(response.context["active_count"], 2)
		self.assertEqual(response.context["inactive_count"], 2)
		self.assertEqual(response.context["availability"], 50)
	
	def test9(self):
		'''track_usage'''
		l1 = Library.objects.get(name="lib1")
		l = LibraryPlate.objects.get(barcode="xyz2", library=l1)
		login = self.client.login(username='staff', password='password')
		response = self.client.get(reverse('track-usage', args=[l.id, '2021-02-03', "general-view"]))
		self.assertEqual(response.status_code, 200)
		self.assertEqual(len(response.context["change_dates"]), 3)
		self.assertEqual(len(response.context["compounds"]), 4)
		self.assertEqual(response.context["active_count"], 3)
		self.assertEqual(response.context["inactive_count"], 1)
		self.assertEqual(response.context["opened"], 1)
		self.assertEqual(response.context["columns"], ['0' + str(i) for i in range(1, 10)] + [str(i) for i in range(10, 25)])
		self.assertEqual(response.context["rows"], [char for char in string.ascii_uppercase][0:16])
		self.assertEqual(response.context["main_id"], "general-view")
		self.assertEqual(response.context["switch_view"], "graphic view")

#	def test10(self):
#		'''proposals'''
#		l1 = Library.objects.get(name="lib1")
#		l = LibraryPlate.objects.get(barcode="xyz2", library=l1)
#		login = self.client.login(username='staff', password='password')
#		response = self.client.get(reverse('proposal'))
#		self.assertEqual(response.status_code, 200)
