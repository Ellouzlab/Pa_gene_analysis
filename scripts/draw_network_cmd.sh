# Factors
nodes=25
min_genomes=2
id=0.8
cov=0.8
# all virulence
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/vfdb_pa_entero.csv' \
--color family  \
--aggregate species \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/species/vfdb' 

# shell
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/pa_shell_entero.csv' \
--color family  \
--aggregate species \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/species/shell' 

# virulence shell
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/pa_vfdb_shell_entero.csv' \
--color family  \
--aggregate species \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/species/shell_vfdb' 

# cloud
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/pa_cloud_entero.csv' \
--color family  \
--aggregate species \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/species/cloud' 

# virulence cloud
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/vfdb_cloud_pa_entero.csv' \
--color family  \
--aggregate species \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/species/cloud_vfdb' 



#####Aggregating by Genus
#virulence
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/vfdb_pa_entero.csv' \
--color family  \
--aggregate genus \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--agglomerans \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/genus/vfdb' 

# shell
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/pa_shell_entero.csv' \
--color family  \
--aggregate genus \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--agglomerans \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/genus/shell' 

# virulence shell
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/pa_vfdb_shell_entero.csv' \
--color family  \
--aggregate genus \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--agglomerans \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/genus/shell_vfdb' 

# cloud
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/pa_cloud_entero.csv' \
--color family  \
--aggregate genus \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--agglomerans \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/genus/cloud' 

# virulence cloud
python draw_network_advanced.py \
--taxonomy '/home/sulman/Desktop/ppanggolin_analysis_final/data/entero_figs/entero_db_mod.csv' \
--hits_csv '/home/sulman/Desktop/ppanggolin_analysis_final/data/vfdb_partition/vfdb_cloud_pa_entero.csv' \
--color family  \
--aggregate genus \
--inclusive \
--min_genomes 2 \
--id $id --cov $cov \
--max_nodes $nodes \
--agglomerans \
--output '/home/sulman/Desktop/ppanggolin_analysis_final/data/network_figs/genus/cloud_vfdb' 