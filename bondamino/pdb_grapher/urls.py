from django.urls import path
from . import views

app_name = 'pdb_grapher'

urlpatterns = [
    path('', views.index, name='index'),
    path('language/', views.set_language, name='set_language'),
    path('account/', views.account, name='account'),
    path('account/theme/', views.set_theme, name='set_theme'),
    path('pdb/<int:pk>/', views.detail, name='detail'),
    path('pdb/<int:pk>/models/', views.model_data, name='model_data'),
    path('pdb/<int:pk>/download/', views.download, name='download'),
    path('pdb/<int:pk>/notes/', views.note_create, name='note_create'),
    path('pdb/<int:pk>/notes/<int:note_pk>/edit/', views.note_edit, name='note_edit'),
    path('pdb/<int:pk>/notes/<int:note_pk>/delete/', views.note_delete, name='note_delete'),
    path('pdb/<int:pk>/notes/<int:note_pk>/groups/<int:group_pk>/remove/', views.note_remove_from_group,
         name='note_remove_from_group'),
    path('pdb/<int:pk>/edit/', views.edit, name='edit'),
    path('pdb/<int:pk>/delete/', views.delete, name='delete'),
    path('pdb/<int:pk>/groups/<int:group_pk>/remove/', views.remove_from_group, name='remove_from_group'),
]
