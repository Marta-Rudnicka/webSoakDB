from django.urls import path
from . import views

urlpatterns = [
    
    path('selection/home/', views.index),
    path('selection/selection/', views.index),
    path('selection/summary/', views.index),
    path('compounds/<str:type>/<int:pk>/<int:proposal_id>/', views.index2),
    path('selection/', views.index),
]
