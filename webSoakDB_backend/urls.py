from rest_framework import routers
from .api import ExampleViewSet
from . import views

from django.contrib import admin
from django.urls import include, path

router = routers.DefaultRouter()
router.register('api/example', ExampleViewSet, 'example')

#urlpatterns = router.urls

urlpatterns = [
    path('', include(router.urls)),
    #path('data_test/', views.data_test, name="data_test"),
    #path('uploads/library_upload_form/library_upload_errors/', views.library_upload_errors, name="library_upload_errors"),
    path('uploads/library_upload_form/', views.upload_user_library, name="library_upload_form"),
    path('uploads/subset_upload_form/', views.upload_user_subset, name="subset_upload_form"),
    path('dummy/', views.dummy, name="dummy"),
]
