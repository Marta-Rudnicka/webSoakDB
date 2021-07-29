from django.shortcuts import render
from django.views.decorators.csrf import ensure_csrf_cookie
from django.contrib.auth.decorators import login_required

@login_required
@ensure_csrf_cookie
def index(request):
    return render(request, 'frontend/index.html')

def index2(request, type, pk, proposal_id):
    return render(request, 'frontend/index.html')
