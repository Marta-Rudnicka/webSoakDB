from rest_framework import serializers
from .models import IspybAuthorization, Library, LibraryPlate, Compounds, Protein, CrystalPlate, Project, LibrarySubset

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

class CompoundSerializer(serializers.Serializer):
	class Meta:
		model = Compounds
		fields = ['id', 
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
				"tpsa"]
		
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
	description = serializers.CharField()
	subsets = LibrarySubsetSerializer(many=True)
	

class CrystalPlateSerializer(serializers.ModelSerializer):
	class Meta:
		model = CrystalPlate
		fields = ['name', 'drop_volume', 'plate_type', 'crystals']
		depth=2

class IspybAuthorizationSerializer(serializers.ModelSerializer):
	class Meta:
		model = IspybAuthorization
		fields = ['id', 'project', 'users', 'proposal_visit']

class ProjectListSerializer(serializers.Serializer):
	id = serializers.IntegerField()
	#proposal = serializers.CharField(max_length=32)
	#libraries = serializers.PrimaryKeyRelatedField(many=True, read_only=True)
	libraries = LibrarySerializer(many=True)
	subsets = LibrarySubsetSerializer(many=True)	
	auth = IspybAuthorizationSerializer(many=True)



class LibrarySubsetSerializer(serializers.ModelSerializer):
	class Meta:
		model = LibrarySubset
		fields=["id", "name", "library", "origin"]
		depth= 2

#CUSTOM SERIALIZERS FOR RUNNING STATS

class CompoundsStatSerializer(serializers.ModelSerializer):
	class Meta:
		model = Compounds
		fields = ['id', 'code', 'smiles', "log_p", "mol_wt", "heavy_atom_count", "heavy_atom_mol_wt", "nhoh_count", "no_count", "num_h_acceptors", "num_h_donors", "num_het_atoms", "num_rot_bonds", "num_val_electrons", "ring_count", "tpsa"]
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
