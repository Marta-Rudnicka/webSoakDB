from django.contrib import admin

# Register your models here.

from .models import Library, Compounds, SourceWell, LibraryPlate, Protein, Project, LibrarySubset, Preset, SpaCompound, Crystal, CrystalPlate, SWStatuschange, PlateOpening



class SourceWellAdmin(admin.ModelAdmin):
    list_display = ("id", "library_plate", "well", "compound", "concentration", "active", "deactivation_date", "status_changes")

class LibraryPlateAdmin(admin.ModelAdmin):
    list_display = ("id", "library", "barcode", "last_tested", "current")

class LibraryAdmin(admin.ModelAdmin):
    list_display = ("id", "name", "public", "for_industry")
    
class CompoundsAdmin(admin.ModelAdmin):
	 #list_display = ("code", "smiles")
	 list_display = ["id", "code", "smiles", "log_p", "mol_wt", "heavy_atom_count", "heavy_atom_mol_wt", "nhoh_count", "no_count", "num_h_acceptors", "num_h_donors", "num_het_atoms", "num_rot_bonds", "num_val_electrons", "ring_count", "tpsa"]
	 list_per_page = 800

class ProjectAdmin(admin.ModelAdmin):
    #filter_horizontal = ['auth']
    list_display = ("protein", "industry_user")

class CrystalPlateAdmin(admin.ModelAdmin):
    list_display = ("id", "name", "drop_volume", "plate_type")

class LibrarySubsetAdmin(admin.ModelAdmin):
	list_display = ("id", "name", "library", "origin")

class PresetAdmin(admin.ModelAdmin):
	list_display = ("id", "name", "description")

class SWStatuschangeAdmin(admin.ModelAdmin):
	list_display = ("id", "source_well", "date", "activation")

class PlateOpeningAdmin(admin.ModelAdmin):
	list_display = ("id", "plate", "date", "reason")

	 
admin.site.register(Compounds, CompoundsAdmin)    
admin.site.register(SourceWell, SourceWellAdmin)
admin.site.register(LibraryPlate, LibraryPlateAdmin)
admin.site.register(Library, LibraryAdmin)
admin.site.register(Project,ProjectAdmin)
admin.site.register(Protein)
admin.site.register(LibrarySubset, LibrarySubsetAdmin)
admin.site.register(Preset, PresetAdmin)
admin.site.register(SpaCompound)
admin.site.register(Crystal)
admin.site.register(CrystalPlate, CrystalPlateAdmin)
admin.site.register(SWStatuschange, SWStatuschangeAdmin)
admin.site.register(PlateOpening, PlateOpeningAdmin)
