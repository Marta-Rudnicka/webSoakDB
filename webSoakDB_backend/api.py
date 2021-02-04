from .models import Example
from rest_framework import viewsets, permissions
from .serializers import ExampleSerializer


# Crystal ViewSet
class ExampleViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = Example.objects.all()
    permissions_classes = [
        permissions.AllowAny
    ]
    serializer_class = ExampleSerializer

