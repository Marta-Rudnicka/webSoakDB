from rest_framework import serializers
from .models import Library, LibraryPlate, Compounds, SourceWell, Protein, Preset, CrystalPlate, Proposals, LibrarySubset

class LibrarySerializer(serializers.ModelSerializer):
	class Meta:
		model = Library
		fields = ['id', 'name', 'for_industry', 'public', 'plates']
		depth = 1

class LibraryPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['id', 'library', 'name', 'current', 'size']
		depth = 1

#class LibraryPlateWithCompoundsSerializer(serializers.ModelSerializer):
#	class Meta:
#		model = LibraryPlate
#		fields = ['library', 'name', 'current', 'compounds', 'size']


class CurrentPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['id', 'library', 'name', 'size']
		depth = 1

class CompoundSerializer(serializers.Serializer):
	class Meta:
		model = Compounds
		fields = ['id', 'code', 'smiles', 'properties']
		
class LibrarySubsetSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	library = LibrarySerializer()
	name = serializers.CharField(max_length = 64)
	origin = serializers.CharField(max_length = 256)
	size = serializers.IntegerField()

'''
class SourceWellSerializer(serializers.ModelSerializer):
	class Meta:
		model = SourceWell
		fields = ['library_plate', 'well', 'compound', 'concentration']
		depth = 1

class SourceWellSerializer(serializers.Serializer):
	well = serializers.CharField(max_length=4)
	compound = CompoundSerializer()
	concentration = serializers.IntegerField()
'''		
		
class ProteinSerializer(serializers.ModelSerializer):
	class Meta:
		model = Protein
		fields = ['id', 'proposal', 'name', 'space_group', 'a', 'b', 'c', 'alpha', 'beta', 'gamma']

 
class PresetSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	name = serializers.CharField(max_length=32)
	description = serializers.CharField(max_length=256)
	subsets = LibrarySubsetSerializer(many=True)

class CrystalPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = CrystalPlate
		fields = ['name', 'drop_volume', 'plate_type', 'crystals']
		depth=2

class ProposalListSerializer(serializers.Serializer):
	name = serializers.CharField(max_length=32)
	#libraries = serializers.PrimaryKeyRelatedField(many=True, read_only=True)
	libraries = LibrarySerializer(many=True)
	subsets = LibrarySubsetSerializer(many=True)	

'''class LibrarySubsetSerializer(serializers.ModelSerializer):
	class Meta:
		model LibrarySubset
		fields'''

#CUSTOM SERIALIZERS FOR RUNNING STATS

class CompoundsStatSerializer(serializers.ModelSerializer):
	class Meta:
		model = Compounds
		fields = ['id', 'code', 'smiles', 'properties']
		depth: 1

class SourceWellStatSerializer(serializers.Serializer):
	compound = CompoundsStatSerializer()
	well = serializers.CharField(max_length=4)
	concentration = serializers.IntegerField()

class LibraryPlateStatSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	name = serializers.CharField(max_length=32)
	library = LibrarySerializer()
	compounds = SourceWellStatSerializer(many=True)

class LibrarySubsetStatSerializer(LibraryPlateStatSerializer):
	compounds = CompoundsStatSerializer(many=True)

#CUSTOM SERIALIZERS FOR COMPOUND SELECTION APP
class ProposalUpdateSerializer(serializers.ModelSerializer):
	class Meta:
		model = Proposals
		fields = ["name", "libraries", "subsets"]
		lookup_field = "name"
