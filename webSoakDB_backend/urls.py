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
    path('data_test/', views.data_test, name="data_test")
]
