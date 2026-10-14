"""Root URL configuration of the Bondamino project.

Project-wide routes are defined here:

    admin/    Django admin (users, work groups, languages, PDB files, notes).
    login/    Sign-in with pdb_grapher.forms.LoginForm (sign-in attempts are throttled there).
    logout/   Sign-out (POST).
    jsi18n/   JavaScript translation catalog for the scripts in pdb_grapher/static.
    ''        Everything else is delegated to pdb_grapher/urls.py (namespace `pdb_grapher`).

New pages of the app belong in pdb_grapher/urls.py, not here. See
https://docs.djangoproject.com/en/5.2/topics/http/urls/ for the URL dispatcher.
"""
from django.contrib import admin
from django.contrib.auth import views as auth_views
from django.urls import path, include
from django.views.i18n import JavaScriptCatalog

from pdb_grapher.forms import LoginForm

# Uploaded PDB files are deliberately not served via MEDIA_URL, but only through the views of
# pdb_grapher, which check the access rights.
urlpatterns = [
    path('admin/', admin.site.urls),
    path('login/', auth_views.LoginView.as_view(template_name='pdb_grapher/login.html', authentication_form=LoginForm,
                                                redirect_authenticated_user=True),
         name='login'),
    path('logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('jsi18n/', JavaScriptCatalog.as_view(packages=['pdb_grapher']), name='javascript-catalog'),
    path('', include('pdb_grapher.urls')),
]
