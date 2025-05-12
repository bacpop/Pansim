# static params
core=1342
pan=4400
freq=0.45 # based on average genome size 2000 genes
ngen=100
npop=1000
mu_core=0.019
core_nuc=$((core * 1000))
pos_lambda=100
neg_lambda=100
rate_genes2=1000
HR_rate=0.0

# changing params - baseline
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}_baseline
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - low prop_positive
prop_positive=0.1
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med prop_positive
prop_positive=0.5
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - high prop_positive
prop_positive=1.0
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - low HGT
prop_positive=-0.1
HGT_rate=0.1
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med HGT
prop_positive=-0.1
HGT_rate=1.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - high HGT
prop_positive=-0.1
HGT_rate=10.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med gene-gain/loss
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=10.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params high gene-gain/loss
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=100.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - low prop2
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.1

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med prop2
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.5

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - high prop2
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=1.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med HGT and med selection
prop_positive=0.5
HGT_rate=1.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med HGT and med gene-gain
prop_positive=0.0
HGT_rate=1.0
rate_genes1=1.0
prop_genes2=0.5

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med selection and med gene-gain
prop_positive=0.5
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.5

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - med selection and med gene-gain and med HGT
prop_positive=0.5
HGT_rate=1.0
rate_genes1=1.0
prop_genes2=0.5

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose

# changing params - strong competition
prop_positive=-0.1
HGT_rate=0.0
rate_genes1=1.0
prop_genes2=0.0

outpref=prop_positive_${prop_positive}_HGT_rate_${HGT_rate}_rate_genes1_${rate_genes1}_prop_genes2_${prop_genes2}_competition_10000
./pansim/target/release/pansim --n_gen ${ngen} --pop_size ${npop} --core_size ${core_nuc} --pan_genes ${pan} --core_genes ${core} --avg_gene_freq ${freq} --threads 4 --core_mu ${mu_core} --HR_rate ${HR_rate} --HGT_rate ${HGT_rate} --rate_genes1 ${rate_genes1} --rate_genes2 ${rate_genes2} --prop_genes2 ${prop_genes2} --prop_positive ${prop_positive} --pos_lambda ${pos_lambda} --neg_lambda ${neg_lambda} --outpref ${outpref} --verbose --competition --competition_strength 10000
