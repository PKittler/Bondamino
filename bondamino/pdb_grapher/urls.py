from django.urls import path
from . import views

app_name = 'pdb_grapher' # Namespace for URLs

urlpatterns = [
    path('', views.index_view, name='index'),
    # No separate ‘upload’ path required, is handled by index_view (POST)
    path('select/<uuid:pdb_file_id>/', views.select_options_view, name='select_options'),
    path('process/<uuid:pdb_file_id>/', views.process_pdb_view, name='process_pdb'),
    path('results/<uuid:graph_id>/', views.show_results_view, name='show_results'),
]