# Generated by Django 3.1.5 on 2021-08-20 15:16

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Example',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('smiles', models.CharField(blank=True, db_index=True, max_length=255, null=True, unique=True)),
            ],
        ),
    ]
