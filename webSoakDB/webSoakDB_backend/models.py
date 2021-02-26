# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from __future__ import unicode_literals

from django.db import models

import os


class Example(models.Model):
    smiles = models.CharField(
        max_length=255, blank=True, null=True, db_index=True, unique=True)
    # Uncomment if using xcdb
    # class Meta:
    #    app_label = 'xchem_db'
    #    db_table = 'compounds'
