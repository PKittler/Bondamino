from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', include('pdb_grapher.urls')),
]

if settings.DEBUG:
    # serve media data only in debug mode
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)