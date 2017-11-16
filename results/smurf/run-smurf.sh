SMURF_DIR="$HOME/Ballmer/SMURF"
LUPIENGROUP_DIR="/mnt/work1/users/lupiengroup/People"
PEAKS="regions_LNCaP-DHS_H3K27ac-catalogue.bed"

sh $SMURF_DIR/smurf.sh \
    -o "$LUPIENGROUP_DIR/hawleyj/Wittenberg/results/smurf" \
    -f "vcf" \
    -v "$LUPIENGROUP_DIR/parisa/PCa_project/mutations_from_BoutrosLab/perSampleVCFs/" \
    -s "$SMURF_DIR/hg19_snp147Common.sorted.bed" \
    -p "$SMURF_DIR/Gencodev24_Promoters_NEW.bed" \
    -r "$PEAKS" \
    -m "regionsamples"
