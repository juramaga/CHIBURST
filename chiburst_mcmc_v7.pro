 ;This program is intended to find the best fit to provided spectral
;and photometric data of star-forming galaxies from a grid of SED models (Groves et al. 2008). Besides,
;it produces chi^2 maps of the parameter space to look for
;degeneracies among parameters and evaluate local minima of the chi^2 distribution.

;CALLING SEQUENCE:
;IDL>  CHIBURST, wave, flux, error, distance, source_no

;PARAMETERS:
;wave        = An array containing the nomimal wavelenths of the observed filters.
;flux        = The measured filter fluxes (in Jy)
;error       = The observational errors in the fluxes (in Jy)
;distance    = The distance to the source, in Mpc
;source_no   = An identifier for the source


;OUTPUT:
;The output of the code is a grid with the chi^2 values por each model
;for the given data. Additionally, the user can get chi^2 maps of the
;desired parameters.


;========================================================================
;COMP_CHI2: This function calculates the chi^2 value for a given
;model and photometry dataset

;CALLING SEQUENCE
;IDL> comp_chi2(fluxes, errors, weight)
; fluxes: The model fluxes at the wavelengths given by 'wave'
; errors: The associated errors
; weight: The weight to be given to each datapoint (uniform, etc.)
;========================================================================

FUNCTION COMP_CHI2, obs_fluxes, obs_errors, model_fluxes, weight

common stats, n_data, log_obs_var2

n_data = n_elements(obs_fluxes) ; number of datapoints

log_obs_flux = alog10(obs_fluxes) - 0.5*(1.0/alog(10.))*(1.0/(obs_fluxes^2.0))*obs_errors*obs_errors ;Unbiased flux in log space

log_obs_var2 = (((1.0/alog(10.))*(1.0/obs_fluxes))^2.0)*obs_errors*obs_errors ;Unbiased observed variances

log_mod_flux = alog10(model_fluxes) ;Model flux in log space

chi_squared = (1.0/n_data)*total(((log_obs_flux - log_mod_flux)*(log_obs_flux - log_mod_flux))/(2.0*weight*(log_obs_var2))) ; reduced chi squared

if chi_squared ne 'Inf' then PDF = exp(-chi_squared) else PDF = 0.0 ;This is the calculated probability

pdf_str = {chi2:chi_squared, pdf:PDF}

return, pdf_str
 
END


;========================================================================
;SED_INTERPOL: This function calculates the interpolation between two
;given gridpoints of the Groves et al. SED models.

;CALLING SEQUENCE
;IDL> 
;========================================================================

FUNCTION SED_INTERPOL, pars
common grid, ALL_Z_MODELS1, CPARAM1, METAL1, F_PDR1, PK1, UCHII_05, UCHII_10, UCHII_20, OLDSTAR, FREQ_HZ, DIFFUSE_EMISSION
UCHII = [[UCHII_05],[UCHII_10],[UCHII_20]]
OLDSTARS = OLDSTAR[2:4,*]
DIFFUSE = DIFFUSE_EMISSION[2:4,6,*]


FOR i=0,n_elements(metal1)-2 DO BEGIN
    
    IF pars[0] GE metal1[i] AND pars[0] LT metal1[i+1] THEN BEGIN
        
        f_met=(pars[0]-metal1[i])/(metal1[i+1]-metal1[i])
        
        ind_met=i+f_met
        
    ENDIF
    
ENDFOR



FOR i=0,n_elements(cparam1)-2 DO BEGIN

    
    IF pars[1] GE cparam1[i] AND pars[1] LT cparam1[i+1] THEN BEGIN
        
        f_c=(pars[1]-cparam1[i])/(cparam1[i+1]-cparam1[i])
        
        ind_cparam=i+f_c
        
    ENDIF
    
ENDFOR



FOR i=0,n_elements(pk1)-2 DO BEGIN
    
    IF pars[2] GE pk1[i] AND pars[2] LT pk1[i+1] THEN BEGIN
        
        f_pk=(pars[2]-pk1[i])/(pk1[i+1]-pk1[i])
        
        ind_pk=i+f_pk
        
    ENDIF
    
ENDFOR


FOR i=0,n_elements(f_pdr1)-2 DO BEGIN
    
    IF pars[3] GE f_pdr1[i] AND pars[3] LT f_pdr1[i+1] THEN BEGIN
        
        f_f_pdr=(pars[3]-f_pdr1[i])/(f_pdr1[i+1]-f_pdr1[i])
        
        ind_f_pdr=i+f_f_pdr
        
    ENDIF

ENDFOR

UCHII_SED = ((1-f_met)*UCHII[*,fix(ind_met)] + f_met*UCHII[*,1+fix(ind_met)])/freq_hz

OLDSTAR_SED = ((1-f_met)*OLDSTARS[fix(ind_met),*] + f_met*OLDSTARS[1+fix(ind_met),*])

DIFFUSE_SED = ((1-f_met)*DIFFUSE[fix(ind_met),0,*] + f_met*DIFFUSE[1+fix(ind_met),0,*])   ;We use the an average ISM field of 10 times G0, appropriate for 'normal' galaxies where the ISM is not extremely heated by old populations.


int_met0A=(1-f_met)*all_z_models1[fix(ind_met),fix(ind_cparam),fix(ind_pk),0,*]+f_met*all_z_models1[1+fix(ind_met),fix(ind_cparam),fix(ind_pk),0,*]

int_met0B=(1-f_met)*all_z_models1[fix(ind_met),1+fix(ind_cparam),fix(ind_pk),0,*]+f_met*all_z_models1[1+fix(ind_met),1+fix(ind_cparam),fix(ind_pk),0,*]

int_met0C=(1-f_met)*all_z_models1[fix(ind_met),fix(ind_cparam),1+fix(ind_pk),0,*]+f_met*all_z_models1[fix(1+ind_met),fix(ind_cparam),1+fix(ind_pk),0,*]

int_met0D=(1-f_met)*all_z_models1[fix(ind_met),1+fix(ind_cparam),1+fix(ind_pk),0,*]+f_met*all_z_models1[fix(1+ind_met),1+fix(ind_cparam),1+fix(ind_pk),0,*]

int_met1A=(1-f_met)*all_z_models1[fix(ind_met),fix(ind_cparam),fix(ind_pk),1,*]+f_met*all_z_models1[1+fix(ind_met),fix(ind_cparam),fix(ind_pk),1,*]

int_met1B=(1-f_met)*all_z_models1[fix(ind_met),1+fix(ind_cparam),fix(ind_pk),1,*]+f_met*all_z_models1[1+fix(ind_met),1+fix(ind_cparam),fix(ind_pk),1,*]

int_met1C=(1-f_met)*all_z_models1[fix(ind_met),fix(ind_cparam),1+fix(ind_pk),1,*]+f_met*all_z_models1[fix(1+ind_met),fix(ind_cparam),1+fix(ind_pk),1,*]

int_met1D=(1-f_met)*all_z_models1[fix(ind_met),1+fix(ind_cparam),1+fix(ind_pk),1,*]+f_met*all_z_models1[fix(1+ind_met),1+fix(ind_cparam),1+fix(ind_pk),1,*]



int_cparam1=(1.0-f_c)*int_met0A+f_c*int_met0B

int_cparam2=(1.0-f_c)*int_met0C+f_c*int_met0D

int_cparam3=(1.0-f_c)*int_met1A+f_c*int_met1B

int_cparam4=(1.0-f_c)*int_met1C+f_c*int_met1D




int_pk1=(1.0-f_pk)*int_cparam1+f_pk*int_cparam2

int_pk2=(1.0-f_pk)*int_cparam3+f_pk*int_cparam4



final_sed=(1-f_f_pdr)*int_pk1+f_f_pdr*int_pk2


structure_seds = {basic:final_sed, uchii:UCHII_SED, oldstars:OLDSTAR_SED, diffuse:DIFFUSE_SED}

RETURN, structure_seds

END




;========================================================================
;CONFIDENCE: Calculates the areas of 68.3% (1-sigma) and
;            90% confidence and translates them into 
;            PDF map contours.
;========================================================================

FUNCTION CONFIDENCE, PDF_MAP

; 1-SIGMA:

conf_contours=0.0*findgen(2)
val_contour=0.0*findgen(n_elements(pdf_map))

FOR k=0,n_elements(pdf_map)-1 DO BEGIN
    indo=array_indices(pdf_map,k)
    val_contour[k]=pdf_map[indo[0],indo[1]]
ENDFOR

index_order=sort(val_contour)
val_contour=val_contour[index_order]
val_contour=reverse(val_contour)

l=0

REPEAT BEGIN
    
    area=0.0
    
    FOR j=0,n_elements(pdf_map)-1 DO BEGIN
        indi=array_indices(pdf_map,j)
        IF (pdf_map[indi[0],indi[1]] GT val_contour[l]) THEN area=area+pdf_map[indi[0],indi[1]]
    ENDFOR
    
    l=l+1
    
ENDREP UNTIL area GT 0.683

conf_contours[0]=val_contour[l-1] ;1-sigma contour

; 90% Confidence:

l=0

REPEAT BEGIN
    
    
    area=0.0
    FOR j=0,n_elements(pdf_map)-1 DO BEGIN
        indi=array_indices(pdf_map,j)
        IF (pdf_map[indi[0],indi[1]] GT val_contour[l]) THEN area=area+pdf_map[indi[0],indi[1]]
    ENDFOR
    
    
    l=l+1
    
ENDREP UNTIL area GT 0.9

conf_contours[1]=val_contour[l-1] ;90% confidence contour


 RETURN,conf_contours
 
END


;========================================================================
;CONFIDENCE_1D: Calculates the areas of 68.3% (1-sigma) and
;               90% confidence and translates them into 
;               error bars.
;========================================================================

FUNCTION CONFIDENCE_1D, HISTOGRAM

; 1-SIGMA:

conf_contours=0.0*findgen(2)
val_contour=histogram


;FOR k=0,n_elements(histogram)-1 DO BEGIN
;    val_contour[k]=histogram[k]
;ENDFOR

index_order=sort(val_contour)
val_contour=val_contour[index_order]
val_contour=reverse(val_contour)

l=0

REPEAT BEGIN
    
    area=0.0
    
    FOR j=0,n_elements(histogram)-1 DO BEGIN
        IF (histogram[j] GT val_contour[l]) THEN area=area+histogram[j]
    ENDFOR
    
    l=l+1
    
ENDREP UNTIL area GT 0.683

conf_contours[0]=val_contour[l-1] ;1-sigma contour

; 90% Confidence:

l=0

REPEAT BEGIN
    
    
    area=0.0
    FOR j=0,n_elements(histogram)-1 DO BEGIN
        IF (histogram[j] GT val_contour[l]) THEN area=area+histogram[j]
    ENDFOR
    
    
    l=l+1
    
ENDREP UNTIL area GT 0.9

conf_contours[1]=val_contour[l-1] ;90% confidence contour


 RETURN,conf_contours
 
END

;========================================================================
;REBIN_DATA: Rebins model SED grid to the resolution of the observatuions. 
;========================================================================

FUNCTION REBIN_DATA,wave_0,wave_ref,flux_0,error_0

ind_order = sort(wave_0)   ;Sort wavelengths in ascending order
wave = wave_0[ind_order]
flux = flux_0[ind_order]
error = error_0[ind_order]


ind_order1 = sort(wave_ref)
wavita = wave_ref[ind_order1]


flux_binned = fltarr(n_elements(wavita))
error_binned = fltarr(n_elements(wavita))

it=where(wavita ge wave[0] and wavita le wave[n_elements(wave)-1])

wavita = wavita[it]

flux_binned = flux_binned[it]

error_binned = error_binned[it]

FOR m=0,n_elements(wavita)-2 DO BEGIN
    
    bin=where(wave ge wavita[m] and wave lt wavita[m+1])
    IF (mean(bin) eq -1.0) THEN BEGIN
        flux_binned[m] = flux_binned[m-1]
        error_binned[m] = error_binned[m-1]
    ENDIF ELSE BEGIN
        bin=bin-n_elements(bin)/2.0
        flux_binned[m]=mean(flux[bin])
        error_binned[m]=sqrt(total((error[bin])^2.0))/n_elements(bin)   ; Propagation of error within the bin
        ;error_binned[m]=mean(error[bin])                                 ; SIMPLE ERROR
    ENDELSE
    
ENDFOR

flux_binned[n_elements(wavita)-1] = flux_binned[n_elements(wavita)-2]
error_binned[n_elements(wavita)-1] = error_binned[n_elements(wavita)-2]

str_binned = {wave_bin:wavita, flux_bin:flux_binned, error_bin:error_binned}

RETURN, str_binned

END




;========================================================================
;REBIN_PHOTOMETRY: Extracts model SED points corresponding to the
;observed photometry 
;========================================================================

FUNCTION REBIN_PHOTOMETRY,wave_0,wave_ref

ind_order = sort(wave_0)   ;Sort wavelengths in ascending order
wave = wave_0[ind_order]
;flux = flux_0[ind_order]
;error = error_0[ind_order]


ind_order1 = sort(wave_ref)
wavita = wave_ref[ind_order1]


;flux_binned = fltarr(n_elements(wavita))
;error_binned = fltarr(n_elements(wavita))

wavota = fltarr(n_elements(wave))



FOR i=0,n_elements(wave)-1 DO BEGIN
    
    it=fix(mean(where(abs(alog10(wavita)-alog10(wave[i])) lt 0.0075)))
    wavota[i] = wavita[it]

ENDFOR

RETURN, wavota

END



;========================================================================
;CHIBURST: Main program
;========================================================================

PRO CHIBURST, wave_phot1, flux_phot_Jy1, error_phot_Jy1, distance, name
common grid, ALL_Z_MODELS1, CPARAM1, METAL1, F_PDR1, PK1, UCHII_05, UCHII_10, UCHII_20, OLDSTAR, FREQ_HZ, DIFFUSE_EMISSION
common stats, n_data, log_obs_var2
common plot, wave_obs_binned,flux_obs_binned,error_obs_binned,wave_model_binned,flux_model_bf_binned,flux_model_binned_init,out_dir,wave_micron,sed4_bf,flux_obs_erg_s_binned,error_obs_erg_s_binned,factor,frequita2,frequita1,frequita3,basico,embebido,viejas,difuso


;-----------------------------------------------------------------------
;DEFINITIONS & CALLS
;-----------------------------------------------------------------------


astrolib

device,decompose=0

loadct,39


zeros = where(flux_phot_Jy1 gt 0.0)  ; Get rid of negative fluxes

wave_phot = wave_phot1[zeros]

flux_phot_Jy = flux_phot_Jy1[zeros]

error_phot_Jy = error_phot_Jy1[zeros]


out_dir = '/data/irac12/jmartine/chiburst/results/elbaz/'

dir_seds = '/data/irac12/jmartine/chiburst/brent_models/' ;Directory where model seds are

restore,'/data/irac12/jmartine/chiburst/frequencies.save' ;For conversion to frequencies

frequita1 = freq_hz

frequita2 = interpol(freq_hz, wave_micron, wave_phot)

distance_cm = distance*3.08567758E24

num_data2 = n_elements(wave_phot)




;-----------------------------------------------------------------------
;MODULE 0:
;This module converts observed fluxes and errors to energy units
;------------------------------------------------------------------------

factor = (1.0e-23)*4.0*!pi*distance_cm
factor = factor*distance_cm

flux_phot = flux_phot_Jy*factor
flux_phot_erg_s = flux_phot*frequita2

error_phot = error_phot_Jy*factor
error_phot_erg_s = error_phot*frequita2

lum_tot = int_tabulated(wave_phot,flux_phot_erg_s,/sort,/double)  ; Aproximate total luminosity of the galaxy, in erg/s

SFR_24 = ((2.04E-43)*mean(flux_phot_erg_s[where(abs(24.0-wave_phot) le 0.1)]))*4.0 ;SFR_24 AS DERIVED FROM CALZETTI LAW (IN M_SUN/YR) (modified by a factor of 4)

print,'SFR24 = ',SFR_24


;THIS CONVERSION IS TO ESTIMATE AN INITIAL GUESS FOR THE OLD POPULATIONS FROM THE RATIO BETWEEN OPTICAL AND IR FLUXES
SFR_OLD_FACT1 = mean(flux_phot_Jy[where(wave_phot gt 1.0 and wave_phot lt 3.0)])   
SFR_OLD_FACT2 = mean(flux_phot_Jy[where(wave_phot gt 75.0 and wave_phot lt 120.0)])
SFR_OLD_FACT3 = mean(flux_phot_Jy[where(wave_phot gt 0.1 and wave_phot lt 0.9)])

print,SFR_OLD_FACT1
print,SFR_OLD_FACT2
print,SFR_OLD_FACT3

FACT_OLD = 300.0*((SFR_OLD_FACT1)^(1./3.))
FACT_AV = 10.0*(SFR_OLD_FACT3/SFR_OLD_FACT1)

print,FACT_OLD, FACT_AV

;Errors that are too small are corrected:

;frac_error = error_phot_Jy/flux_phot_Jy
;help,frac_error
;print,frac_error
;error_small = where(frac_error lt 0.01)


;error_phot_Jy[error_small] = 0.01*flux_phot_Jy[error_small]


;-----------------------------------------------------------------------
;MODULE 1:
;This first reads the model SEDS at the provided
;filters, for a selected aperture. The read values are in erg/s/Hz,
;except UCHII regions, which are in erg/s
;------------------------------------------------------------------------


;Now we read the model fluxes

restore, dir_seds + 'all_z_models.save'  ; Models for all metallicities

all_z_models1 = all_z_models
cparam1 = cparam
metal1 = metal
f_pdr1 = f_pdr
pk1 = pk

restore, dir_seds + 'UCHIIdata.save'          ; UCHII regions

restore, dir_seds + 'Oldstardata2.save'       ; Old Populations

restore, dir_seds + 'Attenuationdata.save'    ; Attenuation data

restore, dir_seds + 'Diffuse_emission.save'        ; Diffuse emission



;-----------------------------------------------------------------------
;MODULE 2: 
;Given the observed spectra and/or photometry and the distance to the source, this
;module interpolates the grid of models to any value inside the
;parameter space, and rebins the data tot he model resolution
;-----------------------------------------------------------------------


;First the random number generator

n = 100000        ; Number of iterations
;x = 0              ; Initial value of candidate


metallicity = fltarr(1)
compactness = fltarr(1)
pressure = fltarr(1)
old_populations = fltarr(1)
sf_rate = fltarr(1)
mass_embedded = fltarr(1)
embedded_fraction = fltarr(1)
pdr_fraction = fltarr(1)
diffuse_em = fltarr(1)
av = fltarr(1)
chisq_list = fltarr(1)
chisq_list[0] = -1
pdf_list = fltarr(1)
pdf_list[0] = -1
aprobacion=fltarr(1)
;sed_out = fltarr(1,1,1,1,1800)
   


; NOW WE DO THE BINNING


    phot_wave_model = rebin_photometry(wave_phot,wave_micron)        ; Model Wavelengths of the photometry

    wave_obs_binned = [phot_wave_model]                              ; Final observed wavelengths, rebinned

    flux_obs_binned = [flux_phot_Jy]                                 ; Final observed fluxes, rebinned

    flux_obs_erg_s_binned = [flux_phot_erg_s]                        ; Final observed fluxes (erg/s), rebinned

    error_obs_binned = [error_phot_Jy]                               ; Final observed errors, rebinned

    error_obs_erg_s_binned = [error_phot_erg_s]                      ; Final observed errors (erg/s), rebinned

    match,wave_obs_binned,wave_micron,suba,subb,epsilon=0.0001       ; Corresponding model bins to the data

    wave_model_binned = wave_micron[subb]                            ; Final modeled wavelengths, rebinned


;INITIALIZE PARAMETERS (**********MAYBE BUILD A FUNCTION, WHAT IS HERE
;MUST BE ALMOST IDENTICAL TO TH LOOP BELOW********)

    init_sfr = 0.5*SFR_24                    ;innov_sfr[x]                   ; This is the SFR in the last 10 Myr (SFR_SB). Range as described above
    init_op  = 0.005*FACT_OLD                       ;init_sfr*innov_op[x]           ; This is the SFR from 100Myr ago until 10Myr ago (SFR_OLD). Varies between 1% and 100% of SFR_SB
    ;init_op  = 50.0*init_sfr
    init_embed  = 0.001*init_sfr            ;*innov_embed[0]     ; This is the SFR in the last million year of massive YSOs (SFR_EMB) still heavilly embeded. Varies same as SFR_OLD
    init_pdr  = 0.5                      ;innov_pdr[x]                  ; This is the fraction of PDR
    init_diff  = 0.0001                    ;innov_diff[x]                ; Diffuse contribution. 
    init_comp  = 5.0                     ;innov_comp[x]                ; Compactness parameter
    init_pk  = 5.0                       ;innov_pk[x]                    ; Ambient pressure
    init_metal = 1.0                     ;innov_metal[x]               ; Metallicity
    init_av = 5.0;2.0/FACT_AV                        ;innov_av[x]                     ; Av, varies between 0 and 10


pars_init = [init_metal,init_comp,init_pk,init_pdr]

sed_first_init = sed_interpol(pars_init)

basic_sed_init = sed_first_init.basic/factor ; Model SED in Jy (HII + PDR)

uchii_sed_init = sed_first_init.uchii/factor ; Model SED in Jy (EMBEDDED)

oldstar_sed_init = sed_first_init.oldstars/factor ; Model SED in Jy  (OLD STARS)

diffuse_sed_init = (1.0E23)*sed_first_init.diffuse/1E6

itol = where(wave_micron ge 0.01 and wave_micron le 1000.)
wavitol = wave_micron[itol]

lum_guess = int_tabulated(wavitol,freq_hz[itol]*factor*diffuse_sed_init[itol],/sort,/double)  ;This part normalizes the cold diffuse emission to a fraction of the total galxy luminosity

print,'LUM_GUESS = ',lum_guess
print,'LUM_TOT = ',lum_tot

diffuse_sed_init = init_diff*(lum_tot/lum_guess)*diffuse_sed_init ; Model diffuse dust emission in Jy  

sed1_init = basic_sed_init + init_embed*uchii_sed_init

sed2_init = sed1_init+init_op*oldstar_sed_init

sed3_init = init_sfr*sed2_init

sed4_init = init_diff*diffuse_sed_init + sed3_init*exp(-init_av*ABS) ; FINAL SED, added extinction

flux_model_binned_init = sed4_init[subb] ; Final modeled fluxes, rebinned



;SETTING UP STEPS

rng=Obj_new('RandomNumberGenerator', initialseed)

;The range in SFR is set according to the 24um calibration of Calzetti
;et al. 2012 (a variation of 1.5 orders of magnitude, logarithmically,
;starting at the value given by the Calzetti law)

;*************THE FIRST 3 (SFR, EMBEDDED, OLD POP, MUST VARY LOGARITHMICALLY)************

innov_sfr = rng -> GetRandomNumbers(n)
;innov_sfr = -0.15*(alog10(init_sfr)) + 0.3*(alog10(init_sfr))*innov_sfr
innov_sfr = -0.15 + 0.3*innov_sfr
innov_sfr = 10.0^innov_sfr

innov_pdr = rng -> GetRandomNumbers(n)
innov_pdr = -0.15+0.3*innov_pdr

innov_diff = rng -> GetRandomNumbers(n)
innov_diff = -0.15+0.3*innov_diff
innov_diff = 10.0^innov_diff

innov_embed = rng -> GetRandomNumbers(n)
innov_embed = -0.15+0.3*innov_embed
innov_embed = 10.0^innov_embed

innov_op = rng -> GetRandomNumbers(n) 
innov_op = -0.15+0.3*innov_op
innov_op = 10.0^innov_op


innov_comp = rng -> GetRandomNumbers(n)
innov_comp = -1.25*0.3+0.3*2.5*innov_comp

innov_pk = rng -> GetRandomNumbers(n)
innov_pk = -2.0*0.3+0.3*4.0*innov_pk

innov_metal = rng -> GetRandomNumbers(n)
innov_metal = -0.8*0.3+0.3*1.6*innov_metal

innov_av = rng -> GetRandomNumbers(n)
innov_av = -0.15+0.3*innov_av         ;From exticntion law in Groves et al. DUST_COL = 40.5 corresponds to Av=10
innov_av = 10.0^innov_av


;forprint,wave_obs_binned,flux_obs_binned

window,0   
plot,wave_obs_binned,flux_obs_binned,/xlog,/ylog,xstyle=1,ystyle=1,xr=[0.1,600.],yr=[1E-5,1E6],psym=2
oplot,wave_obs_binned,flux_obs_binned,psym=4, color=250,thick=3,symsize=2
oploterr,wave_obs_binned,flux_obs_binned,error_obs_binned


;print,metal1,cparam1,pk1,f_pdr1


FOR j=0,n-1 DO BEGIN
     
    

    can_sfr = init_sfr*innov_sfr[j]              ; This is the SFR in the last 10 Myr (SFR_SB). Range as described above
    can_op  = init_op*innov_op[j]                ; This is the SFR from 100Myr ago until 10Myr ago (SFR_OLD). Varies between 1% and 100% of SFR_SB
    can_embed  = init_embed*innov_embed[j]       ; This is the SFR in the last million year of massive YSOs (SFR_EMB) still heavilly embeded. Varies same as SFR_OLD
    can_pdr  = init_pdr + innov_pdr[j]           ; This is the fraction of PDR
    can_diff  = init_diff*innov_diff[j]          ; Diffuse contribution. NOT INCLUDED FOR THE MOMENT
    can_comp  = init_comp + innov_comp[j]        ; Compactness parameter
    can_pk  = init_pk + innov_pk[j]              ; Ambient pressure
    can_metal = init_metal + innov_metal[j]      ; Metallicity
    can_av = init_av*innov_av[j]                 ; Av, varies between 0 and 10

    ;print,'HEYY....',alog10(can_sfr),alog10(can_sfr*can_op),alog10(can_sfr*can_embed)

    ;print,can_metal,can_comp,can_pk,can_sfr,can_av,can_diff,can_op,can_embed

    IF ( (can_metal lt 0.4 or can_metal gt 2.0) or (can_comp lt 4.0 or can_comp gt 6.5) or (can_pk lt 4.0 or can_pk gt 8.0) or (can_sfr lt 0.01) or (can_av lt 0.1) or (can_diff lt 0.00001) or (can_op lt 0.001) or (can_embed lt 0.001*init_sfr)) THEN CONTINUE

    

    pars = [can_metal,can_comp,can_pk,can_pdr]

    
    
    IF ( (pars[0] ge metal1[0] and pars[0] le metal1[n_elements(metal1)-1]) AND (pars[1] ge cparam1[0] and pars[1] le cparam1[n_elements(cparam1)-1]) AND (pars[2] ge pk1[0] and pars[2] le pk1[n_elements(pk1)-1]) AND (pars[3] ge 0.0 and pars[3] le 1.0) AND (can_diff ge 0.0 and can_diff le 0.01) AND (can_sfr*can_embed ge 0.0 AND can_sfr*can_embed le 10.0*init_sfr) AND (can_av ge 0.0)) THEN BEGIN
              
        
        sed_first = sed_interpol(pars)
        
        basic_sed = sed_first.basic/factor ; Model SED in Jy (HII + PDR)
        
        uchii_sed = sed_first.uchii/factor ; Model SED in Jy (EMBEDDED)
        
        oldstar_sed = sed_first.oldstars/factor ; Model SED in Jy  (OLD STARS)

        diffuse_sed = (1.0E23)*(lum_tot/lum_guess)*sed_first.diffuse/1E6   ; Model SED in Jy (DIFFUSE)
        
        sed1 = basic_sed + can_embed*uchii_sed
        
        sed2 = sed1+can_op*oldstar_sed
        
        sed3 = can_sfr*sed2
             
        sed4 = can_diff*diffuse_sed + sed3*exp(-can_av*ABS) ; FINAL SED, added extinction
        
        flux_model_binned = sed4[subb] ; Final modeled fluxes, rebinned        
        
        
        
        
;-----------------------------------------------------------------------
;MODULE 3: 
;Given the observed spectra and/or photometry and the distance to the source, this
;module calculates the chi^2 value for selected  models. It uses
;the function comp_chi2, and a Markov chain to explore the parameer space  
;-----------------------------------------------------------------------
        
        
        
        pdf_structure1 = comp_chi2(flux_obs_binned, error_obs_binned, flux_model_binned, 1.0)
        
        pdf_structure2 = comp_chi2(flux_obs_binned, error_obs_binned, flux_model_binned_init, 1.0)
        
        
        prob1 = pdf_structure1.pdf
        
        prob2 = pdf_structure2.pdf
        
        aprob = min([1.0, prob1/prob2])
        
        
    ENDIF ELSE BEGIN

        prob1 = 0.0

        prob2 = 0.0

        aprob = 0.0

    ENDELSE
    
    
    
    IF aprob EQ 'NaN' or aprob EQ 'Inf' or aprob EQ '-Inf' or prob2 EQ 0.0 then aprob = 0.0
    
    u = rng -> GetRandomNumbers(1)
    
    IF (u LT aprob) THEN BEGIN
        
        init_sfr = can_sfr
        sf_rate = [sf_rate,init_sfr]

        init_op = can_op
        old_populations = [old_populations, init_op]

        init_embed = can_embed
        embedded_fraction = [embedded_fraction, init_embed]

        init_pdr = can_pdr
        pdr_fraction = [pdr_fraction,init_pdr]

        init_diff = can_diff
        diffuse_em = [diffuse_em,init_diff]

        init_comp = can_comp
        compactness = [compactness,init_comp]

        init_pk = can_pk
        pressure = [pressure,init_pk]

        init_metal = can_metal
        metallicity = [metallicity,init_metal]

        init_av = can_av
        av = [av,init_av]

        chisq_list = [chisq_list, pdf_structure1.chi2]
        pdf_list = [pdf_list, pdf_structure1.pdf]
        aprobacion=[aprobacion,aprob]

        ;wait,0.5

        ;forprint,can_sfr,can_op,can_embed,can_pdr,can_diff,can_comp,can_metal,can_av

        oplot,wave_model_binned,flux_model_binned,color=10.0+j*0.5,psym=2
                
        ;sed_out = [[sed_out],[sed4]]

    ENDIF


ENDFOR

;oplot,wave_obs_binned,flux_obs_binned,color=0,thick=3,psym=4
;oplot,wave_model_binned,flux_model_binned_init,color=70,thick=3,linestyle=2
oplot,wave_micron,sed4_init,color=70,thick=3
;oplot,wave_micron,diffuse_sed_init,color=250,thick=3


remove,0,metallicity,sf_rate,old_populations,embedded_fraction,pdr_fraction,diffuse_em,compactness,pressure,av,chisq_list,pdf_list,aprobacion

;sed_out = sed_out[*,1:n_elements(chisq_list)]


help,metallicity,sf_rate,old_populations,embedded_fraction,pdr_fraction,diffuse_em,compactness,pressure,av,chisq_list,pdf_list,aprobacion


;FIND THE BEST FIT

ind_min_chi2 = where(chisq_list eq min(chisq_list)); Index of the minimum chi squared in the array

print, 'The minimum chi_squared is: ', min(chisq_list)
print, 'The maximum chi_squared is: ', max(chisq_list)
print, 'The SFR is: ', sf_rate[ind_min_chi2]
sf_rate_bf = mean(sf_rate[ind_min_chi2])
print, 'The SFR over 100 Myr is: ', old_populations[ind_min_chi2]
old_populations_bf = mean(old_populations[ind_min_chi2])
print, 'The SFR over the last Myr is: ', embedded_fraction[ind_min_chi2]
embedded_fraction_bf = mean(embedded_fraction[ind_min_chi2])
print, 'The PDR fraction is: ', pdr_fraction[ind_min_chi2]
pdr_fraction_bf = mean(pdr_fraction[ind_min_chi2])
print, 'The compactness parameter is: ', compactness[ind_min_chi2]
compactness_bf = mean(compactness[ind_min_chi2])
print, 'The ambient pressure is: ', pressure[ind_min_chi2]
pressure_bf = mean(pressure[ind_min_chi2])
print, 'The metallicity is: ', metallicity[ind_min_chi2]
metallicity_bf = mean(metallicity[ind_min_chi2])
print, 'The Av is: ', av[ind_min_chi2]
av_bf = mean(av[ind_min_chi2])
print, 'The fraction of diffuse emission is: ', diffuse_em[ind_min_chi2]
diffuse_bf = mean(diffuse_em[ind_min_chi2])



; Final model SED

pars_bf = [metallicity_bf,compactness_bf,pressure_bf,pdr_fraction_bf]

sed_first_bf = sed_interpol(pars_bf)

basic_sed_bf = sed_first_bf.basic/factor ; Model SED in Jy (HII + PDR)

basico = exp(-av_bf*ABS)*sf_rate_bf*sed_first_bf.basic/factor

uchii_sed_bf = sed_first_bf.uchii/factor ; Model SED in Jy (EMBEDDED)

oldstar_sed_bf = sed_first_bf.oldstars/factor ; Model SED in Jy  (OLD STARS)

embebido = exp(-av_bf*ABS)*sf_rate_bf*embedded_fraction_bf*uchii_sed_bf

diffuse_sed_bf = (1.0E23)*(lum_tot/lum_guess)*sed_first_bf.diffuse/1E6 ; Model SED in Jy  (DIFFUSE)

difuso = transpose(diffuse_bf*diffuse_sed_bf)

sed1_bf = basic_sed_bf + embedded_fraction_bf*uchii_sed_bf

viejas = exp(-av_bf*ABS)*sf_rate_bf*old_populations_bf*oldstar_sed_bf

sed2_bf = sed1_bf+old_populations_bf*oldstar_sed_bf

sed3_bf = sf_rate_bf*sed2_bf

sed4_bf = difuso + sed3_bf*exp(-av_bf*ABS) ; FINAL SED, added extinction

flux_model_bf_binned = sed4_bf[subb] ; Final modeled fluxes, rebinned

oplot,wave_micron,sed4_bf,color=120,thick=3

;oplot,wave_obs_binned,flux_model_bf_binned,color=250,thick=3

;-----------------------------------------------------------------------
;MODULE 4: 
;This module plots the best fit and the derived PDFs for the parameters
;-----------------------------------------------------------------------

param_array = alog10([[sf_rate],[sf_rate*old_populations],[sf_rate*embedded_fraction],[10.0^pdr_fraction],[10.0^compactness],[10.0^pressure],[10.0^metallicity],[av],[diffuse_em]])
best_fit = alog10([sf_rate_bf,sf_rate_bf*old_populations_bf,sf_rate_bf*embedded_fraction_bf,10.0^pdr_fraction_bf,10.0^compactness_bf,10.0^pressure_bf,10.0^metallicity_bf,av_bf,diffuse_bf])
ax_name = ['log SFR!D10!N (M!D!9n!X!N yr!U-1!N)','log SFR!D100!N (M!D!9n!X!N yr!U-1!N)','log SFR!D1!N (M!D!9n!X!N yr!U-1!N)','f!DPDR!N','log Compactness','log Pressure (K cm!U-3!N)','Metallicity (Z!D!9n!X!N)','log A!DV!N (mag)','log f!Ddiff!N']


param_str = {param:param_array, ax_name:ax_name,best_fit:best_fit}

plot_pdfs, param_str, name



END







