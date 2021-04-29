from rest_framework import routers
from .api import ExampleViewSet
from . import views

from django.contrib import admin
from django.urls import include, path

router = routers.DefaultRouter()
router.register('api/example', ExampleViewSet, 'example')

#urlpatterns = router.urls

urlpatterns = [
    path('', views.redirect_to_login),
    path('', include(router.urls)),
    path('dashboard/', views.dashboard),
    path('uploads/library_upload_form/', views.upload_user_library, name="library_upload_form"),
    path('uploads/subset_upload_form/', views.upload_user_subset, name="subset_upload_form"),
    path('dummy/', views.dummy, name="dummy"),
    path('uploads/formatting_help/', views.formatting, name="formatting_help"),
    path('downloads/plate-map/<int:pk>/', views.download_current_plate_map, name="download_plate_map"),
    path('downloads/export-for-soakdb/', views.export_selection_for_soakdb, name="export_for_soakdb"),
]
