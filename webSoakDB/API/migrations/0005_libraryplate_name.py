# Generated by Django 3.1.5 on 2021-08-08 20:09

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0004_auto_20210808_2002'),
    ]

    operations = [
        migrations.AddField(
            model_name='libraryplate',
            name='name',
            field=models.CharField(blank=True, max_length=100, null=True),
        ),
    ]