# Generated by Django 3.1.5 on 2021-04-12 17:03

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0008_auto_20210412_1703'),
    ]

    operations = [
        migrations.AlterField(
            model_name='library',
            name='name',
            field=models.CharField(blank=True, max_length=100, null=True),
        ),
    ]
