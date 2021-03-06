# Generated by Django 3.1.5 on 2021-04-12 20:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0025_solventbatch'),
    ]

    operations = [
        migrations.AddField(
            model_name='batch',
            name='cryo_transfer_vol',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='batch',
            name='expr_conc',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='batch',
            name='soak_vol',
            field=models.FloatField(blank=True, null=True),
        ),
    ]
