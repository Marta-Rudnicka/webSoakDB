from rest_framework import routers
from .api import ExampleViewSet

from django.contrib import admin
from django.urls import include, path

router = routers.DefaultRouter()
router.register('api/example', ExampleViewSet, 'example')

#urlpatterns = router.urls

urlpatterns = [
    path('', include(router.urls)),
]
