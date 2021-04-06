# Generated by Django 3.1.5 on 2021-03-24 11:41

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0006_auto_20210303_1600'),
    ]

    operations = [
        migrations.AddField(
            model_name='sourcewell',
            name='active',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='sourcewell',
            name='deactivation_date',
            field=models.DateField(blank=True, null=True),
        ),
    ]