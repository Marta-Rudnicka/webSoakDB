# Generated by Django 3.1.5 on 2021-07-26 21:10

from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('API', '0002_auto_20210726_2110'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ispybauthorization',
            name='users',
            field=models.ManyToManyField(blank=True, to=settings.AUTH_USER_MODEL),
        ),
    ]
