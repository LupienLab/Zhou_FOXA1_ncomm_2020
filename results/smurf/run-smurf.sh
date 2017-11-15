SMURF_DIR=$HOME/Documents/SMURF/
PEAKS=regions_LNCaP-DHS_H3K27ac-catalogue.bed

sh $SMURF_DIR/smurf.sh \
    -o . \
    -f vcf
    -v $HOME/Lupiengroup/parisa/PCa_project/mutations_from_BoutrosLab/perSampleVCFs/ \
    -s "c" \
    -p $SMURF_DIR/Gencodev24_Promoters_NEW.bed \
    -r $PEAKS \
    -m regionsamples