from django.shortcuts import render
from django.views.decorators.csrf import ensure_csrf_cookie

# Create your views here.
@ensure_csrf_cookie
def index(request):
    return render(request, 'frontend/index.html')

def plate_lookup(request): #, library, plate):
	return render(request, "frontend/plate_lookup.html")#, {
		#"library": library,
		#"plate": plate,
#	})
