# Generated by Django 3.1.5 on 2021-05-04 14:57

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0036_plateopening'),
    ]

    operations = [
        migrations.AddField(
            model_name='compounds',
            name='mol_image',
            field=models.ImageField(blank=True, null=True, upload_to='images/molecules/'),
        ),
    ]
