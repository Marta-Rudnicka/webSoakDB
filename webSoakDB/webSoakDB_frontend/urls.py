from django.urls import path
from . import views

urlpatterns = [
    
    path('home/', views.index),
    path('selection/', views.index),
    path('summary/', views.index),
    path('compounds/<str:type>/<int:pk>/', views.index),
    path('', views.index),
]
