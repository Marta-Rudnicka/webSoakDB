# Generated by Django 3.1.6 on 2021-02-03 17:38

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('API', '0002_auto_20210203_1643'),
    ]

    operations = [
        migrations.CreateModel(
            name='Crystal',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('crystal_name', models.CharField(db_index=True, max_length=255)),
                ('product', models.CharField(blank=True, max_length=255, null=True)),
                ('status', models.CharField(choices=[('PP', 'preprocessing'), ('PD', 'pandda'), ('RE', 'refinement'), ('CC', 'comp_chem'), ('DP', 'deposition')], default='PP', max_length=2)),
                ('well', models.CharField(blank=True, max_length=4, null=True)),
                ('echo_x', models.IntegerField(blank=True, null=True)),
                ('echo_y', models.IntegerField(blank=True, null=True)),
                ('score', models.IntegerField(blank=True, null=True)),
                ('compound', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='API.compounds')),
                ('crystal_plate', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.PROTECT, to='API.crystalplate')),
            ],
            options={
                'db_table': 'crystal',
            },
        ),
        migrations.CreateModel(
            name='DataProcessing',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('auto_assigned', models.TextField(blank=True, null=True)),
                ('cchalf_high', models.FloatField(blank=True, null=True)),
                ('cchalf_low', models.FloatField(blank=True, null=True)),
                ('cchalf_overall', models.FloatField(blank=True, null=True)),
                ('completeness_high', models.FloatField(blank=True, null=True)),
                ('completeness_low', models.FloatField(blank=True, null=True)),
                ('completeness_overall', models.FloatField(blank=True, null=True)),
                ('dimple_mtz_path', models.TextField(blank=True, null=True)),
                ('dimple_pdb_path', models.TextField(blank=True, null=True)),
                ('dimple_status', models.TextField(blank=True, null=True)),
                ('image_path', models.TextField(blank=True, null=True)),
                ('isig_high', models.FloatField(blank=True, null=True)),
                ('isig_low', models.FloatField(blank=True, null=True)),
                ('isig_overall', models.FloatField(blank=True, null=True)),
                ('lattice', models.TextField(blank=True, null=True)),
                ('log_name', models.TextField(blank=True, null=True)),
                ('logfile_path', models.TextField(blank=True, null=True)),
                ('mtz_name', models.TextField(blank=True, null=True)),
                ('mtz_path', models.TextField(blank=True, null=True)),
                ('multiplicity_high', models.FloatField(blank=True, null=True)),
                ('multiplicity_low', models.FloatField(blank=True, null=True)),
                ('multiplicity_overall', models.FloatField(blank=True, null=True)),
                ('original_directory', models.TextField(blank=True, null=True)),
                ('point_group', models.TextField(blank=True, null=True)),
                ('program', models.TextField(blank=True, null=True)),
                ('r_cryst', models.FloatField(blank=True, null=True)),
                ('r_free', models.FloatField(blank=True, null=True)),
                ('r_merge_high', models.FloatField(blank=True, null=True)),
                ('r_merge_low', models.FloatField(blank=True, null=True)),
                ('r_merge_overall', models.FloatField(blank=True, null=True)),
                ('res_high', models.FloatField(blank=True, null=True)),
                ('res_high_15_sigma', models.FloatField(blank=True, null=True)),
                ('res_high_outer_shell', models.FloatField(blank=True, null=True)),
                ('res_low', models.FloatField(blank=True, null=True)),
                ('res_low_inner_shell', models.FloatField(blank=True, null=True)),
                ('res_overall', models.TextField(blank=True, null=True)),
                ('score', models.FloatField(blank=True, null=True)),
                ('spacegroup', models.TextField(blank=True, null=True)),
                ('status', models.TextField(blank=True, null=True)),
                ('unique_ref_overall', models.IntegerField(blank=True, null=True)),
                ('unit_cell', models.TextField(blank=True, null=True)),
                ('unit_cell_vol', models.FloatField(blank=True, null=True)),
                ('crystal_name', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
            ],
            options={
                'db_table': 'data_processing',
            },
        ),
        migrations.CreateModel(
            name='FragalysisLigand',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ligand_name', models.CharField(max_length=255)),
                ('crystallographic_bound', models.FileField(max_length=500, upload_to='')),
                ('lig_mol_file', models.FileField(max_length=500, upload_to='')),
                ('apo_pdb', models.FileField(max_length=500, upload_to='')),
                ('bound_pdb', models.FileField(max_length=500, upload_to='')),
                ('smiles_file', models.FileField(max_length=500, upload_to='')),
                ('desolvated_pdb', models.FileField(max_length=500, upload_to='')),
                ('solvated_pdb', models.FileField(max_length=500, upload_to='')),
                ('pandda_event', models.FileField(blank=True, max_length=500, upload_to='')),
                ('two_fofc', models.FileField(blank=True, max_length=500, upload_to='')),
                ('fofc', models.FileField(blank=True, max_length=500, upload_to='')),
                ('modification_date', models.BigIntegerField()),
            ],
            options={
                'db_table': 'FragalysisLigand',
            },
        ),
        migrations.CreateModel(
            name='Ligand',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.compounds')),
                ('crystal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
                ('fragalysis_ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.fragalysisligand')),
            ],
            options={
                'db_table': 'ligand',
            },
        ),
        migrations.CreateModel(
            name='MiscFiles',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file', models.FileField(max_length=500, upload_to='')),
                ('description', models.TextField()),
            ],
            options={
                'db_table': 'MiscFiles',
            },
        ),
        migrations.CreateModel(
            name='PanddaAnalysis',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pandda_dir', models.CharField(max_length=255, unique=True)),
            ],
            options={
                'db_table': 'pandda_analysis',
            },
        ),
        migrations.CreateModel(
            name='PanddaEvent',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('event', models.IntegerField(blank=True, db_index=True, null=True)),
                ('event_centroid_x', models.FloatField(blank=True, null=True)),
                ('event_centroid_y', models.FloatField(blank=True, null=True)),
                ('event_centroid_z', models.FloatField(blank=True, null=True)),
                ('event_dist_from_site_centroid', models.TextField(blank=True, null=True)),
                ('lig_centroid_x', models.FloatField(blank=True, null=True)),
                ('lig_centroid_y', models.FloatField(blank=True, null=True)),
                ('lig_centroid_z', models.FloatField(blank=True, null=True)),
                ('lig_dist_event', models.FloatField(blank=True, null=True)),
                ('lig_id', models.TextField(blank=True, null=True)),
                ('pandda_event_map_native', models.TextField(blank=True, null=True)),
                ('pandda_event_map_cut', models.TextField(blank=True, null=True)),
                ('pandda_model_pdb', models.TextField(blank=True, null=True)),
                ('pandda_input_mtz', models.TextField(blank=True, null=True)),
                ('pandda_input_pdb', models.TextField(blank=True, null=True)),
                ('ligand_confidence_inspect', models.TextField(blank=True, null=True)),
                ('ligand_confidence', models.TextField(blank=True, null=True)),
                ('comment', models.TextField(blank=True, null=True)),
                ('interesting', models.BooleanField()),
                ('event_status', models.TextField(blank=True, null=True)),
                ('created_date', models.DateTimeField(auto_now_add=True, null=True)),
                ('modified_date', models.DateTimeField(auto_now=True, null=True)),
                ('ligand_confidence_source', models.CharField(choices=[('NA', 'none'), ('SD', 'soak_db'), ('FS', 'fragspect')], default='NA', max_length=2)),
                ('crystal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
                ('data_proc', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.dataprocessing')),
            ],
            options={
                'db_table': 'pandda_event',
            },
        ),
        migrations.CreateModel(
            name='PanddaRun',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('input_dir', models.TextField(blank=True, null=True)),
                ('pandda_log', models.CharField(max_length=255, unique=True)),
                ('pandda_version', models.TextField(blank=True, null=True)),
                ('sites_file', models.TextField(blank=True, null=True)),
                ('events_file', models.TextField(blank=True, null=True)),
                ('pandda_analysis', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.panddaanalysis')),
            ],
            options={
                'db_table': 'pandda_run',
            },
        ),
        migrations.CreateModel(
            name='Reference',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('reference_pdb', models.CharField(default='not_assigned', max_length=255, null=True, unique=True)),
            ],
            options={
                'db_table': 'reference',
            },
        ),
        migrations.CreateModel(
            name='Target',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('target_name', models.CharField(db_index=True, max_length=255, unique=True)),
            ],
            options={
                'db_table': 'target',
            },
        ),
        migrations.CreateModel(
            name='Tasks',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('task_name', models.CharField(db_index=True, max_length=255)),
                ('uuid', models.CharField(db_index=True, max_length=37, unique=True)),
            ],
            options={
                'db_table': 'tasks',
                'unique_together': {('task_name', 'uuid')},
            },
        ),
        migrations.CreateModel(
            name='SoakdbFiles',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('filename', models.CharField(max_length=255, unique=True)),
                ('modification_date', models.BigIntegerField()),
                ('visit', models.TextField()),
                ('status', models.IntegerField(blank=True, null=True)),
                ('proposal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.proposals')),
            ],
            options={
                'db_table': 'soakdb_files',
            },
        ),
        migrations.CreateModel(
            name='ReviewResponses2',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fedid', models.TextField()),
                ('decision_int', models.IntegerField()),
                ('decision_str', models.TextField()),
                ('reason', models.TextField()),
                ('time_submitted', models.IntegerField()),
                ('Ligand_name', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.ligand')),
                ('crystal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
            ],
            options={
                'db_table': 'review_responses_new',
            },
        ),
        migrations.CreateModel(
            name='ReviewResponses',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('fedid', models.TextField()),
                ('decision_int', models.IntegerField()),
                ('decision_str', models.TextField()),
                ('reason', models.TextField()),
                ('time_submitted', models.IntegerField()),
                ('crystal', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
            ],
            options={
                'db_table': 'review_responses',
            },
        ),
        migrations.CreateModel(
            name='Refinement',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('bound_conf', models.CharField(blank=True, max_length=255, null=True, unique=True)),
                ('cif', models.TextField(blank=True, null=True)),
                ('cif_prog', models.TextField(blank=True, null=True)),
                ('cif_status', models.TextField(blank=True, null=True)),
                ('lig_bound_conf', models.TextField(blank=True, null=True)),
                ('lig_cc', models.TextField(blank=True, null=True)),
                ('lig_confidence', models.TextField(blank=True, null=True)),
                ('lig_confidence_int', models.IntegerField(blank=True, null=True)),
                ('lig_confidence_string', models.TextField(blank=True, null=True)),
                ('matrix_weight', models.TextField(blank=True, null=True)),
                ('molprobity_score', models.FloatField(blank=True, null=True)),
                ('mtz_free', models.TextField(blank=True, null=True)),
                ('mtz_latest', models.TextField(blank=True, null=True)),
                ('outcome', models.IntegerField(blank=True, null=True)),
                ('pdb_latest', models.TextField(blank=True, null=True)),
                ('r_free', models.FloatField(blank=True, null=True)),
                ('ramachandran_favoured', models.TextField(blank=True, null=True)),
                ('ramachandran_outliers', models.TextField(blank=True, null=True)),
                ('rcryst', models.FloatField(blank=True, null=True)),
                ('refinement_path', models.TextField(blank=True, null=True)),
                ('res', models.FloatField(blank=True, null=True)),
                ('rmsd_angles', models.TextField(blank=True, null=True)),
                ('rmsd_bonds', models.TextField(blank=True, null=True)),
                ('spacegroup', models.TextField(blank=True, null=True)),
                ('status', models.TextField(blank=True, null=True)),
                ('crystal_name', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
            ],
            options={
                'db_table': 'refinement',
            },
        ),
        migrations.CreateModel(
            name='PanddaSite',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('site', models.IntegerField(blank=True, db_index=True, null=True)),
                ('site_aligned_centroid_x', models.FloatField(blank=True, null=True)),
                ('site_aligned_centroid_y', models.FloatField(blank=True, null=True)),
                ('site_aligned_centroid_z', models.FloatField(blank=True, null=True)),
                ('site_native_centroid_x', models.FloatField(blank=True, null=True)),
                ('site_native_centroid_y', models.FloatField(blank=True, null=True)),
                ('site_native_centroid_z', models.FloatField(blank=True, null=True)),
                ('pandda_run', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.panddarun')),
            ],
            options={
                'db_table': 'pandda_site',
                'unique_together': {('pandda_run', 'site')},
            },
        ),
        migrations.CreateModel(
            name='PanddaEventStats',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('one_minus_bdc', models.FloatField(blank=True, null=True)),
                ('cluster_size', models.IntegerField(blank=True, null=True)),
                ('glob_corr_av_map', models.FloatField(blank=True, null=True)),
                ('glob_corr_mean_map', models.FloatField(blank=True, null=True)),
                ('loc_corr_av_map', models.FloatField(blank=True, null=True)),
                ('loc_corr_mean_map', models.FloatField(blank=True, null=True)),
                ('z_mean', models.FloatField(blank=True, null=True)),
                ('z_peak', models.FloatField(blank=True, null=True)),
                ('b_factor_scaled', models.FloatField(blank=True, null=True)),
                ('high_res', models.FloatField(blank=True, null=True)),
                ('low_res', models.FloatField(blank=True, null=True)),
                ('r_free', models.FloatField(blank=True, null=True)),
                ('r_work', models.FloatField(blank=True, null=True)),
                ('ref_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_b', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_ln_dev', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_ln_dev_z', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_ln_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_ln_rmsd_z', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_below_four_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_below_four_rmsd_z', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_above_four_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_above_four_rmsd_z', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_rmsd_all', models.FloatField(blank=True, null=True)),
                ('wilson_scaled_rmsd_all_z', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_ln_dev', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_ln_dev_z', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_ln_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_ln_rmsd_z', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_below_four_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_below_four_rmsd_z', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_above_four_rmsd', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_above_four_rmsd_z', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_rmsd_all', models.FloatField(blank=True, null=True)),
                ('wilson_unscaled_rmsd_all_z', models.FloatField(blank=True, null=True)),
                ('resolution', models.FloatField(blank=True, null=True)),
                ('map_uncertainty', models.FloatField(blank=True, null=True)),
                ('obs_map_mean', models.FloatField(blank=True, null=True)),
                ('obs_map_rms', models.FloatField(blank=True, null=True)),
                ('z_map_kurt', models.FloatField(blank=True, null=True)),
                ('z_map_mean', models.FloatField(blank=True, null=True)),
                ('z_map_skew', models.FloatField(blank=True, null=True)),
                ('z_map_std', models.FloatField(blank=True, null=True)),
                ('scl_map_mean', models.FloatField(blank=True, null=True)),
                ('scl_map_rms', models.FloatField(blank=True, null=True)),
                ('event', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.panddaevent')),
            ],
            options={
                'db_table': 'pandda_event_stats',
            },
        ),
        migrations.AddField(
            model_name='panddaevent',
            name='pandda_run',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.panddarun'),
        ),
        migrations.AddField(
            model_name='panddaevent',
            name='refinement',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.refinement'),
        ),
        migrations.AddField(
            model_name='panddaevent',
            name='site',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.panddasite'),
        ),
        migrations.CreateModel(
            name='MetaData',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('Site_Label', models.CharField(max_length=255)),
                ('new_smiles', models.TextField(blank=True)),
                ('alternate_name', models.CharField(blank=True, max_length=255)),
                ('pdb_id', models.CharField(blank=True, max_length=255)),
                ('fragalysis_name', models.CharField(max_length=255, unique=True)),
                ('original_name', models.CharField(max_length=255)),
                ('Ligand_name', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.fragalysisligand')),
            ],
            options={
                'db_table': 'metadata',
            },
        ),
        migrations.AddField(
            model_name='ligand',
            name='target',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.target'),
        ),
        migrations.CreateModel(
            name='FragalysisTarget',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('open', models.BooleanField()),
                ('target', models.CharField(max_length=255)),
                ('metadata_file', models.FileField(blank=True, max_length=500, upload_to='')),
                ('input_root', models.TextField()),
                ('staging_root', models.TextField()),
                ('biomol', models.FileField(blank=True, max_length=500, upload_to='')),
                ('additional_files', models.ManyToManyField(to='API.MiscFiles')),
            ],
            options={
                'db_table': 'FragalysisTarget',
            },
        ),
        migrations.AddField(
            model_name='fragalysisligand',
            name='fragalysis_target',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.fragalysistarget'),
        ),
        migrations.AddField(
            model_name='crystal',
            name='target',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.target'),
        ),
        migrations.AddField(
            model_name='crystal',
            name='visit',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.soakdbfiles'),
        ),
        migrations.CreateModel(
            name='BadAtoms',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('atomid', models.IntegerField()),
                ('comment', models.TextField()),
                ('Ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.ligand')),
                ('Review', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.reviewresponses2')),
            ],
            options={
                'db_table': 'bad_atoms',
            },
        ),
        migrations.CreateModel(
            name='PanddaStatisticalMap',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('resolution_from', models.FloatField(blank=True, null=True)),
                ('resolution_to', models.FloatField(blank=True, null=True)),
                ('dataset_list', models.TextField()),
                ('pandda_run', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='API.panddarun')),
            ],
            options={
                'db_table': 'pandda_statistical_map',
                'unique_together': {('resolution_from', 'resolution_to', 'pandda_run')},
            },
        ),
        migrations.AlterUniqueTogether(
            name='panddaevent',
            unique_together={('site', 'event', 'crystal', 'pandda_run')},
        ),
#        migrations.CreateModel(
#            name='Dimple',
#            fields=[
#                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
#                ('mtz_path', models.CharField(blank=True, max_length=255, null=True)),
#                ('pdb_path', models.CharField(blank=True, max_length=255, null=True)),
#                ('r_free', models.FloatField(blank=True, null=True)),
#                ('res_high', models.FloatField(blank=True, null=True)),
#                ('status', models.TextField(blank=True, null=True)),
#                ('crystal_name', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='API.crystal')),
#                ('reference', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='API.reference')),
#            ],
#            options={
#                'db_table': 'dimple',
#                'unique_together': {('pdb_path', 'mtz_path')},
#            },
#        ),
        migrations.AlterUniqueTogether(
            name='crystal',
            unique_together={('crystal_name', 'visit', 'compound', 'product')},
        ),
    ]
