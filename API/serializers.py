from rest_framework import serializers
from .models import Library, LibraryPlate, Compounds, SourceWell, Protein, Preset, CrystalPlate, Proposals

class LibrarySerializer(serializers.ModelSerializer):
	class Meta:
		model = Library
		fields = ['id', 'name', 'for_industry', 'public', 'plates']

class LibraryPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['id', 'library', 'name', 'current', 'size', 'compounds']
		depth = 2

class LibraryPlateWithCompoundsSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['library', 'name', 'current', 'compounds', 'size']


class CurrentPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['id', 'library', 'name', 'size']
		depth = 1

class CompoundSerializer(serializers.ModelSerializer):
	class Meta:
		model = Compounds
		fields = ['code', 'smiles', 'molecular_weight'] #add when RDKit is added
'''
class SourceWellSerializer(serializers.ModelSerializer):
	class Meta:
		model = SourceWell
		fields = ['library_plate', 'well', 'compound', 'concentration']
		depth = 1
'''
class SourceWellSerializer(serializers.Serializer):
	well = serializers.CharField(max_length=4)
	compound = CompoundSerializer()
	concentration = serializers.IntegerField()
		
		
class ProteinSerializer(serializers.ModelSerializer):
	class Meta:
		model = Protein
		fields = ['id', 'proposal', 'name', 'space_group', 'a', 'b', 'c', 'alpha', 'beta', 'gamma']

 
class PresetSerializer(serializers.ModelSerializer):
	class Meta:
		model = Preset
		fields = ['name', 'description', 'subsets']
		depth=2

class CrystalPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = CrystalPlate
		fields = ['name', 'drop_volume', 'plate_type', 'crystals']
		depth=2

class ProposalListSerializer(serializers.ModelSerializer):
	class Meta:
		model = Proposals
		fields = ["name", "libraries", "subsets"]
		depth = 2
		lookup_field = "name"

#CUSTOM SERIALIZERS FOR RUNNING STATS

class CompoundsStatSerializer(serializers.Serializer):
	molecular_weight = serializers.FloatField()

class SourceWellStatSerializer(serializers.Serializer):
	compound = CompoundsStatSerializer()
	concentration = serializers.IntegerField()

class LibraryPlateStatSerializer(serializers.Serializer):
	name = serializers.CharField(max_length=4)
	library = LibrarySerializer()
	compounds = SourceWellStatSerializer(many=True)

#CUSTOM SERIALIZERS FOR COMPOUND SELECTION APP
class ProposalUpdateSerializer(serializers.ModelSerializer):
	class Meta:
		model = Proposals
		fields = ["name", "libraries", "subsets"]
		lookup_field = "name"
