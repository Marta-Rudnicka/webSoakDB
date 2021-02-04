from django.shortcuts import render

# Create your views here.
def index(request):
    return render(request, 'frontend/index.html')

def plate_lookup(request): #, library, plate):
	return render(request, "frontend/plate_lookup.html")#, {
		#"library": library,
		#"plate": plate,
#	})
