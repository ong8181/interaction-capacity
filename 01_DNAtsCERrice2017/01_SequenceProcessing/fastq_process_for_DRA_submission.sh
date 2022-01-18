####
#### Preparation for DRA submission of Claident-processed fastq
####

#--------------- For 515F-806R ---------------#
# Rename fastq
#(for f in *ProSTD40000-*; do mv $f ${f/ProSTD40000-*_S/Prokaryote_S}; done)
for f in *RMR-076__*; do mv $f ${f/RMR-076__/RMR-076_}; done
for f in *__515F*; do mv $f ${f/__515F/_515F}; done

# Decompress fastq.gz
gzip -d *.gz

# Rename fastq sequence ID
find ./ -name "*.fastq" | xargs sed -i 's/__RMR-076__/ RMR-076_/g'
find ./ -name "*.fastq" | xargs sed -i 's/__515F/_515F/g'

# Re-compress fastq
gzip *.fastq


#--------------- For ITS ---------------#
# Rename fastq
for f in *RMR-078__*; do mv $f ${f/RMR-078__/RMR-078_}; done
for f in *__ITS1_F_KYO1*; do mv $f ${f/__ITS1_F_KYO1/_ITS1_F_KYO1}; done

# Decompress fastq.gz
gzip -d *.gz

# Rename fastq sequence ID
find ./ -name "*.fastq" | xargs sed -i 's/__RMR-078__/ RMR-078_/g'
find ./ -name "*.fastq" | xargs sed -i 's/__ITS1_F_KYO1/_ITS1_F_KYO1/g'

# Re-compress fastq
gzip *.fastq


#--------------- For COI ---------------#
# Rename fastq
for f in *RMR-099__*; do mv $f ${f/RMR-099__/RMR-099_}; done
for f in *__mlCOIintF*; do mv $f ${f/__mlCOIintF/_mlCOIintF}; done

# Decompress fastq.gz
gzip -d *.gz

# Rename fastq sequence ID
find ./ -name "*.fastq" | xargs sed -i 's/__RMR-099__/ RMR-099_/g'
find ./ -name "*.fastq" | xargs sed -i 's/__mlCOIintF/_mlCOIintF/g'

# Re-compress fastq
gzip *.fastq




#--------------- For Euk1391f-EukBr ---------------#
# Rename fastq
for f in *CMR-002__*; do mv $f ${f/CMR-002__/CMR-002_}; done
for f in *__Euk1391f*; do mv $f ${f/__Euk1391f/_Euk1391f}; done

# Decompress fastq.gz
gzip -d *.gz

# Rename fastq sequence ID
find ./ -name "*.fastq" | xargs sed -i 's/__CMR-002__/ CMR-002_/g'
find ./ -name "*.fastq" | xargs sed -i 's/__Euk1391f/_Euk1391f/g'

# Re-compress fastq
gzip *.fastq

