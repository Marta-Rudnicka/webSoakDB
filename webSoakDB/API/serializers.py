from django.db.models.fields import BooleanField
from rest_framework import serializers
from .models import IspybAuthorization, Library, LibraryPlate, Compounds, Protein, CrystalPlate, Project, LibrarySubset, SourceWell


class LibrarySerializer(serializers.ModelSerializer):
	class Meta:
		model = Library
		fields = ['id', 'name', 'for_industry', 'public', 'plates']
		depth = 1

class LibraryPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['id', 'library', 'barcode', 'current', 'size']
		depth = 1

class CurrentPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibraryPlate
		fields = ['id', 'library', 'barcode', 'size']
		depth = 1

class LibrarySubsetSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	library = LibrarySerializer()
	name = serializers.CharField(max_length = 64)
	origin = serializers.CharField(max_length = 256)
	size = serializers.IntegerField()

class PresetSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	name = serializers.CharField(max_length=32)
	description = serializers.CharField()
	subsets = LibrarySubsetSerializer(many=True)

class IspybAuthorizationSerializer(serializers.ModelSerializer):
	class Meta:
		model = IspybAuthorization
		fields = ['id', 'project', 'users', 'proposal_visit']

class ProjectListSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	libraries = LibrarySerializer(many=True)
	subsets = LibrarySubsetSerializer(many=True)	
	auth = IspybAuthorizationSerializer(many=True)

class CompoundsStatSerializer(serializers.ModelSerializer):
	class Meta:
		model = Compounds
		fields = [
			'id', 
			'code', 
			'smiles', 
			"log_p", 
			"mol_wt", 
			"heavy_atom_count", 
			"heavy_atom_mol_wt", 
			"nhoh_count", 
			"no_count", 
			"num_h_acceptors", 
			"num_h_donors", 
			"num_het_atoms", 
			"num_rot_bonds", 
			"num_val_electrons", 
			"ring_count", 
			"tpsa"
			]
		depth: 1

class SourceWellStatSerializer(serializers.Serializer):
	compound = CompoundsStatSerializer()
	well = serializers.CharField(max_length=4)
	concentration = serializers.IntegerField()


class LibraryPlateStatSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	barcode = serializers.CharField(max_length=32)
	library = LibrarySerializer()
	compounds = SourceWellStatSerializer(many=True)

class LibrarySubsetStatSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	name = serializers.CharField(max_length=32)
	library = LibrarySerializer()
	compounds = CompoundsStatSerializer(many=True)


#CUSTOM SERIALIZERS FOR COMPOUND SELECTION APP
class ProjectUpdateSerializer(serializers.ModelSerializer):
	class Meta:
		model = Project
		fields = ["proposal", "libraries", "subsets"]
		lookup_field = "proposal"


#SERIALIZERS FOR SECRET DATA

class LibraryInProjectSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	name = serializers.CharField(max_length=64)
	public = BooleanField()
	plates = LibraryPlateStatSerializer(many=True)

class ProjectCompoundsSerializer(serializers.Serializer):
	libraries = LibraryInProjectSerializer(many=True)
