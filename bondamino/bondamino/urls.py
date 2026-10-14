"""
URL configuration for bondamino project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.2/topics/http/urls/
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
from django.contrib.auth import views as auth_views
from django.urls import path, include
from django.views.i18n import JavaScriptCatalog

from pdb_grapher.forms import LoginForm

# Hochgeladene PDB-Dateien werden bewusst nicht über MEDIA_URL ausgeliefert,
# sondern nur über die Views von pdb_grapher, die die Zugriffsrechte prüfen.
urlpatterns = [
    path('admin/', admin.site.urls),
    path('login/', auth_views.LoginView.as_view(template_name='pdb_grapher/login.html', authentication_form=LoginForm,
                                                redirect_authenticated_user=True),
         name='login'),
    path('logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('jsi18n/', JavaScriptCatalog.as_view(packages=['pdb_grapher']), name='javascript-catalog'),
    path('', include('pdb_grapher.urls')),
]
