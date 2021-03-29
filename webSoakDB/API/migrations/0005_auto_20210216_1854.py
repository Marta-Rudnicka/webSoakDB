# Generated by Django 3.1.5 on 2021-02-16 18:54

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0004_auto_20210203_1814'),
    ]

    operations = [
        migrations.AddField(
            model_name='soakdbcompound',
            name='crystal',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.PROTECT, related_name='compounds', to='API.crystal'),
        ),
        migrations.CreateModel(
            name='SolventNotes',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('solvent', models.CharField(choices=[('DMSO', 'DMSO'), ('EG', 'Ethylene Glycol')], default='DMSO', max_length=4)),
                ('solvent_concentration', models.FloatField(blank=True, null=True)),
                ('soak_time', models.DurationField(blank=True, null=True)),
                ('cryo', models.CharField(choices=[('none', 'No cryoprotectant'), ('c1', 'Some string'), ('c2', 'Some other string'), ('c3', 'Yet another string')], default='none', max_length=4)),
                ('cryo_concentration', models.FloatField(blank=True, default=None, null=True)),
                ('comments', models.TextField(blank=True, null=True)),
                ('proposal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.proposals')),
            ],
        ),
    ]