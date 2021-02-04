from django.urls import path
from . import views

urlpatterns = [
    path('', views.index),
   # path('<str:library>/<str:plate>', views.plate_lookup, name="plate_lookup")
    path('test/x', views.plate_lookup, name="plate_lookup")
]
