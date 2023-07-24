
Problematic data sets:

```
- F200W
- jw01837_obs003_visit021_NIRCAM_F200W_i2d -- big problem - (jw01837003021_08201_00001_nrca1_cal)
- jw01837_obs003_visit022_NIRCAM_F200W_i2d -- big problem
- jw01837_obs003_visit023_NIRCAM_F200W_i2d -- big problem
- jw01837_obs003_visit024_NIRCAM_F200W_i2d -- big problem
- jw01837_obs004_visit003_NIRCAM_F200W_i2d -- small claws
- jw01837_obs004_visit004_NIRCAM_F200W_i2d -- small claws
- jw01837_obs004_visit005_NIRCAM_F200W_i2d -- small claws
- jw01837_obs004_visit021_NIRCAM_F200W_i2d -- weird background of one detector part

- F277W
- jw01837_obs004_visit001_NIRCAM_F277W_i2d -- small claws
- jw01837_obs004_visit005_NIRCAM_F277W_i2d -- small claws
- jw01837_obs004_visit014_NIRCAM_F277W_i2d -- some claws at bottom right
-- No big issue

- F356W
- jw01837_obs003_visit021_NIRCAM_F356W_i2d -- big problem - (jw01837003021_04201_00001_nrcalong_cal)
- jw01837_obs003_visit022_NIRCAM_F356W_i2d -- big problem
- jw01837_obs003_visit023_NIRCAM_F356W_i2d -- big problem
- jw01837_obs003_visit024_NIRCAM_F356W_i2d -- big problem
- jw01837_obs004_visit001_NIRCAM_F356W_i2d -- too many saturations but probably fine

- F356W
- jw01837_obs003_visit021_NIRCAM_F444W_i2d -- big problem - (jw01837003021_08201_00001_nrcalong_cal)
- jw01837_obs003_visit022_NIRCAM_F444W_i2d -- big problem
- jw01837_obs003_visit023_NIRCAM_F444W_i2d -- big problem
- jw01837_obs003_visit024_NIRCAM_F444W_i2d -- big problem


- ls -1 jw01837003021_08201_0000*_nrc*/calibrated2_cals/*_cal.fits | wc -l # 20
- ls -1 jw01837003022_08201_0000*_nrc*/calibrated2_cals/*_cal.fits | wc -l # 20
- ls -1 jw01837003023_08201_0000*_nrc*/calibrated2_cals/*_cal.fits | wc -l # 20
- ls -1 jw01837003024_08201_0000*_nrc*/calibrated2_cals/*_cal.fits | wc -l # 20
- ls -1 jw01837003024_04201_0000*_nrc*long/calibrated2_cals/*_cal.fits | wc -l # 4


ds9 jw01837003021_08201_0000*_nrc*/calibrated2_cals/*_cal.fits

# manually draw ploygon regions...

```

```
cp -i jw01837003021_08201_00001_nrca3_bad_area.pix.reg      jw01837003022_08201_00001_nrca3_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrca3_bad_area.pix.reg      jw01837003022_08201_00002_nrca3_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrca4_bad_area.pix.reg      jw01837003022_08201_00001_nrca4_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrca4_bad_area.pix.reg      jw01837003022_08201_00002_nrca4_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcb3_bad_area.pix.reg      jw01837003022_08201_00001_nrcb3_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcb3_bad_area.pix.reg      jw01837003022_08201_00002_nrcb3_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcb4_bad_area.pix.reg      jw01837003022_08201_00001_nrcb4_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcb4_bad_area.pix.reg      jw01837003022_08201_00002_nrcb4_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcalong_bad_area.pix.reg   jw01837003022_08201_00001_nrcalong_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcalong_bad_area.pix.reg   jw01837003022_08201_00002_nrcalong_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcblong_bad_area.pix.reg   jw01837003022_08201_00001_nrcblong_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcblong_bad_area.pix.reg   jw01837003022_08201_00002_nrcblong_bad_area.pix.reg


ds9 -lock frame image -scale zscale \
    jw01837003022_08201_00001_nrca1/calibrated2_cals/jw01837003022_08201_00001_nrca1_cal.fits \
    jw01837003022_08201_00002_nrca1/calibrated2_cals/jw01837003022_08201_00002_nrca1_cal.fits \
    jw01837003022_08201_00001_nrca2/calibrated2_cals/jw01837003022_08201_00001_nrca2_cal.fits \
    jw01837003022_08201_00002_nrca2/calibrated2_cals/jw01837003022_08201_00002_nrca2_cal.fits \
    jw01837003022_08201_00001_nrca3/calibrated2_cals/jw01837003022_08201_00001_nrca3_cal.fits \
        -regions load ../jw01837003022_08201_00001_nrca3_bad_area.pix.reg \
    jw01837003022_08201_00002_nrca3/calibrated2_cals/jw01837003022_08201_00002_nrca3_cal.fits \
        -regions load ../jw01837003022_08201_00002_nrca3_bad_area.pix.reg \
    jw01837003022_08201_00001_nrca4/calibrated2_cals/jw01837003022_08201_00001_nrca4_cal.fits \
        -regions load ../jw01837003022_08201_00001_nrca4_bad_area.pix.reg \
    jw01837003022_08201_00002_nrca4/calibrated2_cals/jw01837003022_08201_00002_nrca4_cal.fits \
        -regions load ../jw01837003022_08201_00002_nrca4_bad_area.pix.reg \
    jw01837003022_08201_00001_nrcb1/calibrated2_cals/jw01837003022_08201_00001_nrcb1_cal.fits \
    jw01837003022_08201_00002_nrcb1/calibrated2_cals/jw01837003022_08201_00002_nrcb1_cal.fits \
    jw01837003022_08201_00001_nrcb2/calibrated2_cals/jw01837003022_08201_00001_nrcb2_cal.fits \
    jw01837003022_08201_00002_nrcb2/calibrated2_cals/jw01837003022_08201_00002_nrcb2_cal.fits \
    jw01837003022_08201_00001_nrcb3/calibrated2_cals/jw01837003022_08201_00001_nrcb3_cal.fits \
        -regions load ../jw01837003022_08201_00001_nrcb3_bad_area.pix.reg \
    jw01837003022_08201_00002_nrcb3/calibrated2_cals/jw01837003022_08201_00002_nrcb3_cal.fits \
        -regions load ../jw01837003022_08201_00002_nrcb3_bad_area.pix.reg \
    jw01837003022_08201_00001_nrcb4/calibrated2_cals/jw01837003022_08201_00001_nrcb4_cal.fits \
        -regions load ../jw01837003022_08201_00001_nrcb4_bad_area.pix.reg \
    jw01837003022_08201_00002_nrcb4/calibrated2_cals/jw01837003022_08201_00002_nrcb4_cal.fits \
        -regions load ../jw01837003022_08201_00002_nrcb4_bad_area.pix.reg \
    jw01837003022_08201_00001_nrcalong/calibrated2_cals/jw01837003022_08201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003022_08201_00001_nrcalong_bad_area.pix.reg \
    jw01837003022_08201_00002_nrcalong/calibrated2_cals/jw01837003022_08201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003022_08201_00002_nrcalong_bad_area.pix.reg \
    jw01837003022_08201_00001_nrcblong/calibrated2_cals/jw01837003022_08201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003022_08201_00001_nrcblong_bad_area.pix.reg \
    jw01837003022_08201_00002_nrcblong/calibrated2_cals/jw01837003022_08201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003022_08201_00002_nrcblong_bad_area.pix.reg \

```

```
cp -i jw01837003021_08201_00001_nrca3_bad_area.pix.reg      jw01837003023_08201_00001_nrca3_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrca3_bad_area.pix.reg      jw01837003023_08201_00002_nrca3_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrca4_bad_area.pix.reg      jw01837003023_08201_00001_nrca4_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrca4_bad_area.pix.reg      jw01837003023_08201_00002_nrca4_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcb3_bad_area.pix.reg      jw01837003023_08201_00001_nrcb3_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcb3_bad_area.pix.reg      jw01837003023_08201_00002_nrcb3_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcb4_bad_area.pix.reg      jw01837003023_08201_00001_nrcb4_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcb4_bad_area.pix.reg      jw01837003023_08201_00002_nrcb4_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcalong_bad_area.pix.reg   jw01837003023_08201_00001_nrcalong_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcalong_bad_area.pix.reg   jw01837003023_08201_00002_nrcalong_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcblong_bad_area.pix.reg   jw01837003023_08201_00001_nrcblong_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcblong_bad_area.pix.reg   jw01837003023_08201_00002_nrcblong_bad_area.pix.reg

ds9 -lock frame image -scale zscale \
    jw01837003023_08201_00001_nrca1/calibrated2_cals/jw01837003023_08201_00001_nrca1_cal.fits \
    jw01837003023_08201_00002_nrca1/calibrated2_cals/jw01837003023_08201_00002_nrca1_cal.fits \
    jw01837003023_08201_00001_nrca2/calibrated2_cals/jw01837003023_08201_00001_nrca2_cal.fits \
    jw01837003023_08201_00002_nrca2/calibrated2_cals/jw01837003023_08201_00002_nrca2_cal.fits \
    jw01837003023_08201_00001_nrca3/calibrated2_cals/jw01837003023_08201_00001_nrca3_cal.fits \
        -regions load ../jw01837003023_08201_00001_nrca3_bad_area.pix.reg \
    jw01837003023_08201_00002_nrca3/calibrated2_cals/jw01837003023_08201_00002_nrca3_cal.fits \
        -regions load ../jw01837003023_08201_00002_nrca3_bad_area.pix.reg \
    jw01837003023_08201_00001_nrca4/calibrated2_cals/jw01837003023_08201_00001_nrca4_cal.fits \
        -regions load ../jw01837003023_08201_00001_nrca4_bad_area.pix.reg \
    jw01837003023_08201_00002_nrca4/calibrated2_cals/jw01837003023_08201_00002_nrca4_cal.fits \
        -regions load ../jw01837003023_08201_00002_nrca4_bad_area.pix.reg \
    jw01837003023_08201_00001_nrcb1/calibrated2_cals/jw01837003023_08201_00001_nrcb1_cal.fits \
    jw01837003023_08201_00002_nrcb1/calibrated2_cals/jw01837003023_08201_00002_nrcb1_cal.fits \
    jw01837003023_08201_00001_nrcb2/calibrated2_cals/jw01837003023_08201_00001_nrcb2_cal.fits \
    jw01837003023_08201_00002_nrcb2/calibrated2_cals/jw01837003023_08201_00002_nrcb2_cal.fits \
    jw01837003023_08201_00001_nrcb3/calibrated2_cals/jw01837003023_08201_00001_nrcb3_cal.fits \
        -regions load ../jw01837003023_08201_00001_nrcb3_bad_area.pix.reg \
    jw01837003023_08201_00002_nrcb3/calibrated2_cals/jw01837003023_08201_00002_nrcb3_cal.fits \
        -regions load ../jw01837003023_08201_00002_nrcb3_bad_area.pix.reg \
    jw01837003023_08201_00001_nrcb4/calibrated2_cals/jw01837003023_08201_00001_nrcb4_cal.fits \
        -regions load ../jw01837003023_08201_00001_nrcb4_bad_area.pix.reg \
    jw01837003023_08201_00002_nrcb4/calibrated2_cals/jw01837003023_08201_00002_nrcb4_cal.fits \
        -regions load ../jw01837003023_08201_00002_nrcb4_bad_area.pix.reg \
    jw01837003023_08201_00001_nrcalong/calibrated2_cals/jw01837003023_08201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003023_08201_00001_nrcalong_bad_area.pix.reg \
    jw01837003023_08201_00002_nrcalong/calibrated2_cals/jw01837003023_08201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003023_08201_00002_nrcalong_bad_area.pix.reg \
    jw01837003023_08201_00001_nrcblong/calibrated2_cals/jw01837003023_08201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003023_08201_00001_nrcblong_bad_area.pix.reg \
    jw01837003023_08201_00002_nrcblong/calibrated2_cals/jw01837003023_08201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003023_08201_00002_nrcblong_bad_area.pix.reg \


```


```
cp -i jw01837003021_08201_00001_nrca3_bad_area.pix.reg      jw01837003024_08201_00001_nrca3_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrca3_bad_area.pix.reg      jw01837003024_08201_00002_nrca3_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrca4_bad_area.pix.reg      jw01837003024_08201_00001_nrca4_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrca4_bad_area.pix.reg      jw01837003024_08201_00002_nrca4_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcb3_bad_area.pix.reg      jw01837003024_08201_00001_nrcb3_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcb3_bad_area.pix.reg      jw01837003024_08201_00002_nrcb3_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcb4_bad_area.pix.reg      jw01837003024_08201_00001_nrcb4_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcb4_bad_area.pix.reg      jw01837003024_08201_00002_nrcb4_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcalong_bad_area.pix.reg   jw01837003024_08201_00001_nrcalong_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcalong_bad_area.pix.reg   jw01837003024_08201_00002_nrcalong_bad_area.pix.reg
cp -i jw01837003021_08201_00001_nrcblong_bad_area.pix.reg   jw01837003024_08201_00001_nrcblong_bad_area.pix.reg
cp -i jw01837003021_08201_00002_nrcblong_bad_area.pix.reg   jw01837003024_08201_00002_nrcblong_bad_area.pix.reg


ds9 -lock frame image -scale zscale \
    jw01837003024_08201_00001_nrca1/calibrated2_cals/jw01837003024_08201_00001_nrca1_cal.fits \
    jw01837003024_08201_00002_nrca1/calibrated2_cals/jw01837003024_08201_00002_nrca1_cal.fits \
    jw01837003024_08201_00001_nrca2/calibrated2_cals/jw01837003024_08201_00001_nrca2_cal.fits \
    jw01837003024_08201_00002_nrca2/calibrated2_cals/jw01837003024_08201_00002_nrca2_cal.fits \
    jw01837003024_08201_00001_nrca3/calibrated2_cals/jw01837003024_08201_00001_nrca3_cal.fits \
        -regions load ../jw01837003024_08201_00001_nrca3_bad_area.pix.reg \
    jw01837003024_08201_00002_nrca3/calibrated2_cals/jw01837003024_08201_00002_nrca3_cal.fits \
        -regions load ../jw01837003024_08201_00002_nrca3_bad_area.pix.reg \
    jw01837003024_08201_00001_nrca4/calibrated2_cals/jw01837003024_08201_00001_nrca4_cal.fits \
        -regions load ../jw01837003024_08201_00001_nrca4_bad_area.pix.reg \
    jw01837003024_08201_00002_nrca4/calibrated2_cals/jw01837003024_08201_00002_nrca4_cal.fits \
        -regions load ../jw01837003024_08201_00002_nrca4_bad_area.pix.reg \
    jw01837003024_08201_00001_nrcb1/calibrated2_cals/jw01837003024_08201_00001_nrcb1_cal.fits \
    jw01837003024_08201_00002_nrcb1/calibrated2_cals/jw01837003024_08201_00002_nrcb1_cal.fits \
    jw01837003024_08201_00001_nrcb2/calibrated2_cals/jw01837003024_08201_00001_nrcb2_cal.fits \
    jw01837003024_08201_00002_nrcb2/calibrated2_cals/jw01837003024_08201_00002_nrcb2_cal.fits \
    jw01837003024_08201_00001_nrcb3/calibrated2_cals/jw01837003024_08201_00001_nrcb3_cal.fits \
        -regions load ../jw01837003024_08201_00001_nrcb3_bad_area.pix.reg \
    jw01837003024_08201_00002_nrcb3/calibrated2_cals/jw01837003024_08201_00002_nrcb3_cal.fits \
        -regions load ../jw01837003024_08201_00002_nrcb3_bad_area.pix.reg \
    jw01837003024_08201_00001_nrcb4/calibrated2_cals/jw01837003024_08201_00001_nrcb4_cal.fits \
        -regions load ../jw01837003024_08201_00001_nrcb4_bad_area.pix.reg \
    jw01837003024_08201_00002_nrcb4/calibrated2_cals/jw01837003024_08201_00002_nrcb4_cal.fits \
        -regions load ../jw01837003024_08201_00002_nrcb4_bad_area.pix.reg \
    jw01837003024_08201_00001_nrcalong/calibrated2_cals/jw01837003024_08201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003024_08201_00001_nrcalong_bad_area.pix.reg \
    jw01837003024_08201_00002_nrcalong/calibrated2_cals/jw01837003024_08201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003024_08201_00002_nrcalong_bad_area.pix.reg \
    jw01837003024_08201_00001_nrcblong/calibrated2_cals/jw01837003024_08201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003024_08201_00001_nrcblong_bad_area.pix.reg \
    jw01837003024_08201_00002_nrcblong/calibrated2_cals/jw01837003024_08201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003024_08201_00002_nrcblong_bad_area.pix.reg \

```


```
cp -i jw01837003024_04201_00001_nrcalong_bad_area.pix.reg jw01837003021_04201_00001_nrcalong_bad_area.pix.reg
cp -i jw01837003024_04201_00002_nrcalong_bad_area.pix.reg jw01837003021_04201_00002_nrcalong_bad_area.pix.reg
cp -i jw01837003024_04201_00001_nrcblong_bad_area.pix.reg jw01837003021_04201_00001_nrcblong_bad_area.pix.reg
cp -i jw01837003024_04201_00002_nrcblong_bad_area.pix.reg jw01837003021_04201_00002_nrcblong_bad_area.pix.reg

cp -i jw01837003024_04201_00001_nrcalong_bad_area.pix.reg jw01837003022_04201_00001_nrcalong_bad_area.pix.reg
cp -i jw01837003024_04201_00002_nrcalong_bad_area.pix.reg jw01837003022_04201_00002_nrcalong_bad_area.pix.reg
cp -i jw01837003024_04201_00001_nrcblong_bad_area.pix.reg jw01837003022_04201_00001_nrcblong_bad_area.pix.reg
cp -i jw01837003024_04201_00002_nrcblong_bad_area.pix.reg jw01837003022_04201_00002_nrcblong_bad_area.pix.reg

cp -i jw01837003024_04201_00001_nrcalong_bad_area.pix.reg jw01837003023_04201_00001_nrcalong_bad_area.pix.reg
cp -i jw01837003024_04201_00002_nrcalong_bad_area.pix.reg jw01837003023_04201_00002_nrcalong_bad_area.pix.reg
cp -i jw01837003024_04201_00001_nrcblong_bad_area.pix.reg jw01837003023_04201_00001_nrcblong_bad_area.pix.reg
cp -i jw01837003024_04201_00002_nrcblong_bad_area.pix.reg jw01837003023_04201_00002_nrcblong_bad_area.pix.reg


ds9 -lock frame image -scale zscale \
    jw01837003021_04201_00001_nrcalong/calibrated2_cals/jw01837003021_04201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003021_04201_00001_nrcalong_bad_area.pix.reg \
    jw01837003021_04201_00002_nrcalong/calibrated2_cals/jw01837003021_04201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003021_04201_00002_nrcalong_bad_area.pix.reg \
    jw01837003021_04201_00001_nrcblong/calibrated2_cals/jw01837003021_04201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003021_04201_00001_nrcblong_bad_area.pix.reg \
    jw01837003021_04201_00002_nrcblong/calibrated2_cals/jw01837003021_04201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003021_04201_00002_nrcblong_bad_area.pix.reg \
    \
    jw01837003022_04201_00001_nrcalong/calibrated2_cals/jw01837003022_04201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003022_04201_00001_nrcalong_bad_area.pix.reg \
    jw01837003022_04201_00002_nrcalong/calibrated2_cals/jw01837003022_04201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003022_04201_00002_nrcalong_bad_area.pix.reg \
    jw01837003022_04201_00001_nrcblong/calibrated2_cals/jw01837003022_04201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003022_04201_00001_nrcblong_bad_area.pix.reg \
    jw01837003022_04201_00002_nrcblong/calibrated2_cals/jw01837003022_04201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003022_04201_00002_nrcblong_bad_area.pix.reg \
    \
    jw01837003023_04201_00001_nrcalong/calibrated2_cals/jw01837003023_04201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003023_04201_00001_nrcalong_bad_area.pix.reg \
    jw01837003023_04201_00002_nrcalong/calibrated2_cals/jw01837003023_04201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003023_04201_00002_nrcalong_bad_area.pix.reg \
    jw01837003023_04201_00001_nrcblong/calibrated2_cals/jw01837003023_04201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003023_04201_00001_nrcblong_bad_area.pix.reg \
    jw01837003023_04201_00002_nrcblong/calibrated2_cals/jw01837003023_04201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003023_04201_00002_nrcblong_bad_area.pix.reg \
    \
    jw01837003024_04201_00001_nrcalong/calibrated2_cals/jw01837003024_04201_00001_nrcalong_cal.fits \
        -regions load ../jw01837003024_04201_00001_nrcalong_bad_area.pix.reg \
    jw01837003024_04201_00002_nrcalong/calibrated2_cals/jw01837003024_04201_00002_nrcalong_cal.fits \
        -regions load ../jw01837003024_04201_00002_nrcalong_bad_area.pix.reg \
    jw01837003024_04201_00001_nrcblong/calibrated2_cals/jw01837003024_04201_00001_nrcblong_cal.fits \
        -regions load ../jw01837003024_04201_00001_nrcblong_bad_area.pix.reg \
    jw01837003024_04201_00002_nrcblong/calibrated2_cals/jw01837003024_04201_00002_nrcblong_cal.fits \
        -regions load ../jw01837003024_04201_00002_nrcblong_bad_area.pix.reg \
```












