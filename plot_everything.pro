PRO PLOT_PDFs,param_str,source_name

common plot

!P.THICK=4
!X.THICK=4
!Y.THICK=4
!Z.THICK=4
!P.CHARTHICK=4
!P.CHARSIZE=1.5
!X.CHARSIZE=1.
!Y.CHARSIZE=1.


nbins=50 ;Numbers of bins for the PDF (histogram)

wave_micron1 = wave_micron*(1+redshift)

set_plot,'ps'

device,/landscape,/color,bits=8,filename=out_dir+'best_fit_'+strtrim(string(source_name),2)+'.ps',xsize=35,ysize=16

;First plot the best fit

multiplot
!p.position = [0.07,0.3,0.45,0.95]

;plot,wave_obs_binned,flux_obs_binned,/xlog,/ylog,xstyle=1,ystyle=1,xr=[0.2,600.],yr=[1E-5,1E3],psym=2,xtitle='Wavelength (!7l!X!Nm)',ytitle='Flux (Jy)'

plot,wave_obs_binned,alog10(factor*frequita2*flux_obs_binned/3.839E33),/xlog,xstyle=1,ystyle=1,xr=[0.07,2000.],yr=[8.0,14.5],psym=2,ytitle='log (!7k!X!NF!D!7k!X!N) [L!D!9n!X!N]',xticklayout=0,xtickname=REPLICATE(' ',7)

;oplot,wave_obs_binned,flux_obs_binned,psym=4, color=250,thick=3,symsize=2

;oploterr,wave_obs_binned,alog10(flux_obs_erg_s_binned/3.839E33),alog10(error_obs_erg_s_binned/3.839E33)
;oplot,wave_model_binned,flux_model_binned_init,color=70,thick=3,linestyle=2

;oplot,wave_micron,alog10(factor*frequita1*sed4_bf/3.839E33),color=70,thick=5
;oplot,wave_micron,alog10(factor*frequita2*flux_model_bf_binned/3.839E33),color=70,thick=5

xyouts,0.2,12.0,source_name

it = where(wave_micron1 gt 0.1 and wave_micron1 lt 2000.)

; Print ascii file with best fit SED

wave_print = reverse(wave_micron1[it])
flux_print = reverse(alog10(factor*frequita1[it]*sed4_bf[it]/3.839E33))

openw,1,out_dir+'best_fit_'+strtrim(string(source_name),2)+'.dat'

printf,1,'Wavelength','log lambda*F_lambda'
printf,1,'[microns]','[L_sun]'

help, factor, frequita1, it, wave_print, flux_print

for i=0,n_elements(wave_print)-1 do printf,1,wave_print[i],flux_print[i]

printf,1,''
printf,1,''


oplot,wave_micron1,alog10(factor*frequita1*basico/3.839E33),color=90,thick=4,linestyle=2
oplot,wave_micron1,alog10(factor*frequita1*viejas/3.839E33),color=230,thick=4,linestyle=3
oplot,wave_micron1,alog10(factor*frequita1*embebido/3.839E33),color=210,thick=4,linestyle=4
oplot,wave_micron1,alog10(factor*frequita1*viejisimas/3.839E33),color=40,thick=4,linestyle=2


oplot,wave_obs_binned,alog10(flux_obs_erg_s_binned/3.839E33),psym=4,thick=3,symsize=1.5


oplot,wave_micron1,alog10(factor*frequita1*sed4_bf/3.839E33),color=70,thick=5
;oplot,wave_obs_binned,alog10(factor*frequita2*flux_model_bf_binned/3.839E33),color=70,thick=5

!p.position = [0.07,0.11,0.45,0.3]

help,frequita2,flux_model_bf_binned

plot,wave_obs_binned,alog10(flux_obs_erg_s_binned/3.839E33)-alog10(factor*frequita2*flux_model_bf_binned/3.839E33),xtitle='Wavelength [!7l!X!Nm]',psym=4,xstyle=1,ystyle=1,/xlog,yr=[-1.0,1.0],ytickinterval=1.0,xr=[0.07,2000.]


oplot,[0.0001,10000.],[0.0,0.0],linestyle=2



;Now plot the PDFs

FOR i = 0, 8 DO BEGIN


    ;print,param_str.param[i,*]
    
    histo_param = histogram(param_str.param[*,i],nbins=nbins,locations=param_values) ;The histogram with the PDF is produced here.
    histo_param = histo_param/total(histo_param) ;Normalized probability

    
    lev = confidence_1d(histo_param) ;Probability level for which 1-sigma is attained
    
    limits = param_values[sort(abs(histo_param - lev[0]))] ; Limits of the 1sigma region
    
    limits90 = param_values[sort(abs(histo_param - lev[1]))] ; Limits of the 90% confidence region


    !p.position = [0.50+(i mod 3)*0.18,0.1+0.3*(i/3),0.63+(i mod 3)*0.18,0.3+0.3*(i/3)]
    
    plot,param_values,histo_param,psym=10,xstyle=1,xtitle=param_str.ax_name[i], ytitle='Likelihood',charsize=1.0;,xticks=3,xtickformat='(F6.2)'

   
    

    ;FOR 90% CONFIDENCE

    limits0 = limits90[0]
    l = 1
    ;limits1 = limits0
    
    REPEAT BEGIN
        
        limits1 = limits90[l] 
        l = l+1
        
    ENDREP UNTIL abs(limits0-limits1) gt 25.0*(param_values[n_elements(param_values)-1]-param_values[0])/nbins

  

    pdf_peak = where(histo_param eq max(histo_param))
    pos_peak = param_values[pdf_peak[0]] ; Parameter value where the PDF reaches its maximum.
    

    IF (limits0 lt limits1) and (pos_peak gt limits0 and pos_peak lt limits1) THEN BEGIN 
        
        conf90 = where(param_values ge limits0 and param_values le limits1)
        poly_x = [limits0-0.0001,param_values[conf90],limits1-0.0001]
        poly_y = [0,histo_param[conf90],0]
 
    ENDIF ELSE IF (limits1 lt limits0) and (pos_peak gt limits1 and pos_peak lt limits0) THEN BEGIN
                  
        conf90 = where(param_values ge limits1 and param_values le limits0)
        poly_x = [limits1-0.0001,param_values[conf90],limits0-0.0001]
        poly_y = [0,histo_param[conf90],0]

    ENDIF ELSE BEGIN

        conf90 = where(param_values gt limits0 and param_values lt max(param_values))
        poly_x = [limits0-0.0001,param_values[conf90],max(param_values)-0.0001]
        poly_y = [0,histo_param[conf90],0]
        
    ENDELSE
    
    forprint,poly_x[0],poly_x[n_elements(poly_x)-1]
    forprint,poly_y[0],poly_y[n_elements(poly_y)-1]
    print,limits0
    print,' '
    

    polyfill,poly_x,poly_y,color=250,/line_fill,linestyle=0,orientation=45,spacing=0.05

    
    ;FOR ONE SIGMA

    limits0 = limits[0]
    l = 1
    
    REPEAT BEGIN
        
        limits1 = limits[l] 
        l = l+1
        
    ENDREP UNTIL abs(limits0-limits1) gt 15.0*(param_values[n_elements(param_values)-1]-param_values[0])/nbins
    
    
    IF limits0 lt limits1 THEN BEGIN 
        
        one_sigma = where(param_values ge limits0 and param_values le limits1)
        poly_x = [limits0-0.0001,param_values[one_sigma],limits1-0.0001]
        poly_y = [0,histo_param[one_sigma],0]
        
        limites0 = limits0
        limites1 =limits1
        
    ENDIF ELSE BEGIN
        
        one_sigma = where(param_values ge limits1 and param_values le limits0)
        poly_x = [limits1-0.0001,param_values[one_sigma],limits0-0.0001]
        poly_y = [0,histo_param[one_sigma],0]
        
        limites0 = limits1
        limites1 =limits0
        
    ENDELSE

    polyfill,poly_x,poly_y,color=70,/line_fill,linestyle=0,orientation=90,spacing=0.05

    ;print,param_str.best_fit[i]

    oplot,[param_str.best_fit[i],param_str.best_fit[i]],[0,1],linestyle=2


    ; Also, print PDF:
   
    printf,1,param_str.ax_name[i]
    printf,1,'Best fit : ',param_str.ax_name[i]+' = '+param_str.best_fit[i]
    printf,1,'Peak of PDF : ',param_str.ax_name[i]+' = '+param_values[pdf_peak[0]]
    for j=0, n_elements(param_values)-1 do printf,1,param_values[j],histo_param[j]
    printf,1,''
    printf,1,''

    
ENDFOR

close,1

device,/close

save,param_str,param_values,histo_param,filename=out_dir+'best_fit_'+strtrim(string(source_name),2)+'.save'

cleanplot
set_plot,'x'


END




