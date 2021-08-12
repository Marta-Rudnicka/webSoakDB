from django.shortcuts import render
from django.core.exceptions import ObjectDoesNotExist
from django.http import HttpResponseRedirect
from django.contrib.admin.views.decorators import staff_member_required
from webSoakDB_stack.settings import MEDIA_ROOT
from . import inv_helpers as h
from .. import forms as forms
from API.models import Project, IspybAuthorization
from django.contrib.auth.models import User

#VIEWS HANDLING GET REQUESTS
@staff_member_required
def index(request):
	return render(request, "inventory/inventory-index.html")

@staff_member_required
def projects(request):
	old_project_form = forms.OldProjectForm()
	new_project_form = forms.NewProjectForm()
	return render(request, "inventory/projects.html", {
		"old_project_form" : old_project_form, 
		"new_project_form" : new_project_form
		})

@staff_member_required
def manage_users(request):
	new_staff_member = forms.AddStaffUserFrom()
	new_power_user = forms.AddPowerUserFrom()

	staff_members = User.objects.filter(is_staff=True)
	power_users = []
	print('staff_members: ', staff_members)
	return render(request, "inventory/users.html", {
		"staff_members": staff_members,
		"power_users" : power_users,
		"new_staff_member" : new_staff_member,
		"new_power_user" : new_power_user
	})

@staff_member_required
def add_staff_member(request):
	form = forms.AddStaffUserFrom(request.POST)
	if request.method == "POST":
		if form.is_valid():
			username = form.cleaned_data['username_staff']
			first_name = form.cleaned_data['first_name_staff']
			last_name = form.cleaned_data['last_name_staff']
			email = form.cleaned_data['email_staff']

			try:
				user = User.objects.get(username=username)
				if first_name:
					user.first_name = first_name
				if last_name:
					user.last_name = last_name
				if email:
					user.email = email
				user.is_staff=True
				user.save()
			except ObjectDoesNotExist:
				user = User.objects.create(
					username=username,
					first_name=first_name,
					last_name = last_name,
					email=email,
					is_staff=True
				)
				user.set_unusable_password()
				user.save()
			return HttpResponseRedirect('../users')
		else:
			return render(request, "inventory/errors.html", {'error_log': [form.errors, form.non_field_errors]})

@staff_member_required
def remove_from_staff(request):
	if request.method == "POST":
		user_id = request.POST.get('user_id')
		try:
			user = User.objects.get(pk=user_id)
			user.is_staff=False
			user.save()
		except ObjectDoesNotExist:
			return render(request, "inventory/errors.html", {'error_log': "Application error: trying to edit non-existent user"})
		return HttpResponseRedirect('../users')

@staff_member_required
def proposal(request):
	form = forms.OldProjectForm(request.POST)
	if request.method == "POST":	
		if form.is_valid():
			try:
				pr = form.cleaned_data['proposal']
				project = h.get_project_by_proposal(pr) #Project.objects.get(proposal=pr)
				subsets = h.get_subsets_with_availability(set(project.subsets.all()))
				visits = IspybAuthorization.objects.filter(proposal_visit__startswith=pr)
								
				return render(request, "inventory/proposal.html", {
					'project' : project,
					'proposal' : pr,
					'subsets': subsets, 
					'visits' : visits
					#'visit_form' : forms.AddVisitForm()
					})
			except(ObjectDoesNotExist):
				return render(request, "inventory/errors.html", {'error_log': ['Proposal not found']})
				
		else:
			return render(request, "inventory/errors.html", {'error_log': [form.errors, form.non_field_errors]})

def dummy(request):
	return render(request, "webSoakDB_backend/dummy.html");

@staff_member_required
def add_project(request):
	if request.method == "POST":
		form = forms.NewProjectForm(request.POST)
		if form.is_valid():
			proposal = request.POST["proposal"]
			title = request.POST["title"]
			if request.POST.get("industry_user", False):
				industry_user = True
			else:
				industry_user = False

			first_visit = proposal + '-1'

			new_project = Project.objects.create(industry_user=industry_user)
			
			new_auth = IspybAuthorization.objects.create(project=title, proposal_visit=first_visit)			
			new_project.auth.add(new_auth)
			new_auth.save()
			
			return HttpResponseRedirect('../projects')
		else:
			return render(request, "inventory/errors.html")	

def add_visit(request):
	if request.method == "POST":
		number = request.POST['number']
		proposal = request.POST['proposal']
		proposal_visit = proposal + '-' + str(number)
		project = Project.objects.filter(auth__proposal_visit__startswith=proposal)[0]
		project_title = project.auth.all()[0].project + '-' + str(number)
		new_auth = IspybAuthorization.objects.create(project=project_title, proposal_visit=proposal_visit)	
		project.auth.add(new_auth)
		new_auth.save()

		visits = IspybAuthorization.objects.filter(proposal_visit__startswith=proposal)
		return render(request, "inventory/proposal.html", {
				'project' : project,
				'proposal' : proposal,
				'subsets': h.get_subsets_with_availability(set(project.subsets.all())), 
				'visits' : visits
				})
