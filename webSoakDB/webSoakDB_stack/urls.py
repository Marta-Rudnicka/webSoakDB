"""webSoakDB_stack URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from django.conf.urls import url
from django.conf import settings
from django.conf.urls.static import static
import django_cas_ng.views

urlpatterns = [
    
    path('admin/', admin.site.urls),
    path('api/', include('API.urls')),
    path('inventory/', include('inventory.urls')),
    #path('accounts/', include('django.contrib.auth.urls')),
    path('', include('webSoakDB_backend.urls')),
    url(r"^accounts/login/", django_cas_ng.views.LoginView.as_view(), name="cas_ng_login"),
    url(r"^accounts/logout/", django_cas_ng.views.LogoutView.as_view(), name="cas_ng_logout"),
    url(
        r"^accounts/callback$",
        django_cas_ng.views.CallbackView.as_view(),
        name="cas_ng_proxy_callback",
    ),
    path('', include('webSoakDB_frontend.urls')),
    
]+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

