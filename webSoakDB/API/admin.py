from django.contrib import admin

# Register your models here.

from .models import Library, Compounds, SourceWell, LibraryPlate, Protein, Proposals, LibrarySubset, Preset, SoakDBCompound, Crystal, CrystalPlate



class SourceWellAdmin(admin.ModelAdmin):
    list_display = ("id", "library_plate", "well", "compound", "concentration", "active", "deactivation_date")

class LibraryPlateAdmin(admin.ModelAdmin):
    list_display = ("id", "library", "name", "last_tested", "current")

class LibraryAdmin(admin.ModelAdmin):
    list_display = ("id", "name", "public", "for_industry")
    
class CompoundsAdmin(admin.ModelAdmin):
	 list_display = ("code", "smiles")
	 list_per_page = 800

class ProposalsAdmin(admin.ModelAdmin):
    list_display = ("name", "protein", "industry_user")

class CrystalPlateAdmin(admin.ModelAdmin):
    list_display = ("id", "name", "drop_volume", "plate_type")

class LibrarySubsetAdmin(admin.ModelAdmin):
	list_display = ("id", "name", "library", "origin")

class PresetAdmin(admin.ModelAdmin):
	list_display = ("id", "name", "description")
	 
admin.site.register(Compounds, CompoundsAdmin)    
admin.site.register(SourceWell, SourceWellAdmin)
admin.site.register(LibraryPlate, LibraryPlateAdmin)
admin.site.register(Library, LibraryAdmin)
admin.site.register(Proposals,ProposalsAdmin)
admin.site.register(Protein)
admin.site.register(LibrarySubset, LibrarySubsetAdmin)
admin.site.register(Preset, PresetAdmin)
admin.site.register(SoakDBCompound)
admin.site.register(Crystal)
admin.site.register(CrystalPlate, CrystalPlateAdmin)
