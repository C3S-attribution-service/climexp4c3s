all climexp.a wrappers.a selectdate lastvalid select_min_years netcdf2ascii cumul catnc diff_hist_nat synthesis setundef count_missing outliers wetbulb_field compute_wetbulb testwetbulb printbigtable ecmwf_times getchance extend_series plotdaily attributefield patchfield diamond2year extremeseries extreme.h makesnow polygon2box quantiles_field quantiles_series averageseries polygon2mask makeleap makeweek grads2nc nc2varlist hurricane_vecchi fieldsignificance fix_undef subtractfield patchseries trendfield fillin del_dimension spectrum subfieldseries timeshift geowind multifit seriesensanomal month2lead flattennc flattennc_dec convertmetadata seriesanomal fieldclim averagefieldspace average_ensemble averagefield_ensemble statmodel extractseries extractfield untransform coordinates2kml list2kml RPS ROCscoremap rocdeb Briar rocmap roc difffield runningmean operate series getunits regionverification verification ar1 txt2dat svd fieldcorrelate daily2longerfield diffdat getval yearly2shorter yearly2shorterfield get_depth dat2grads ctl2dat findmax scaleseries runningmoments maskseries getmomentsfield normdiff get_index get_index_mask netcdf2dat grads2ascii wave sstoi2dat getnperyear daily2longer climatology correlatefield correlatefieldfield filteryearfield filteryearseries filtermonthseries correlate dat2nc selectyear grib2nc plotdat describefield stationlist eof autocor lomb patternfield histogram attribute attribute_amoeba transform maketestfile testprog testerfcc test_fitcross testfitgevcov testshiftseries month2string season2string annual2string gen_time F2f90 install clean:	. $(PVM_ARCH)
	cd $(PVM_ARCH); make $@
