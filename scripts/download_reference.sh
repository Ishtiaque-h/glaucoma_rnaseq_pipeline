# Create a reference folder
OUTDIR="references"
mkdir -p $OUTDIR

# === Ensembl example (GRCm39, release variableized) ===
ENSEMBL_RELEASE=112
BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}

# Transcriptome (cDNA + ncRNA for comprehensive quant)
curl -O ${BASE}/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
curl -O ${BASE}/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz
cat Mus_musculus.GRCm39.cdna.all.fa.gz Mus_musculus.GRCm39.ncrna.fa.gz > "$OUTDIR/mmus_transcripts.fa.gz"
rm -rf Mus_musculus.GRCm39.cdna.all.fa.gz
rm -rf Mus_musculus.GRCm39.ncrna.fa.gz

# Annotation (GTF)
curl -O ${BASE}/gtf/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz
mv Mus_musculus.GRCm39.${ENSEMBL_RELEASE}.gtf.gz "$OUTDIR/mmus.gtf.gz"
