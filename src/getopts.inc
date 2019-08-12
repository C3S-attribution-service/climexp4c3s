        real decor,kscut,minfac,minindx,maxindx,mindata,maxdata,         &
     &      pminindx,pmaxindx,pmindata,pmaxdata,lon1,lon2,lat1,lat2,     &
     &      lev1,lev2,altlon1,altlon2,altlat1,altlat2,altlev1,altlev2,   &
     &      confidenceinterval,restrain,minfacsum,xyear,biasadd,biasmul, &
     &      biasrt
        integer lag1,lag2,m1,m2,yr1,yr2,yr1a,yr2a,yr1b,yr2b,lsum,        &
     &      lsum2,lsel,                                                  &
     &      ndiff,ndiff2,indxuse,avex,avey,nfittime,minnum,intertype,    &
     &      mdiff,mdiff2,nens1,nens2,nyrwindow,lead1,lead2,irunvar,      &
     &      noisetype,debias,normalization,fitfunc,altavex,altavey,      &
     &      ncrossvalidate,day0,add_option
        logical logscale,sqrtscale,lrank,ldetrend,lsubtract,lwrite,      &
     &      lincl(19),lks,dump,plot,lconting,fix2,lnormsd,lbootstrap,    &
     &      anom,composite,ensemble,squarescale,cubescale,               &
     &      twothirdscale,leads,                                         &
     &      logfield,sqrtfield,lfitnoise,lrandom,lrmse,lchangesign,      &
     &      lmakeensfull,lweb,lstandardunits,lallobs,lstartstop,         &
     &      lnomissing,lensanom,lnooverlap,lmae,lincludelast
        character oper*1,strindx(19)*20,plotfile*255,                    &
     &      indexfiles(10:19)*255,plotrunfile*255,corrmonths*8,          &
     &      namestring*80,bbfile*255
        common /c_getopts/ decor,minfac,kscut,                           &
     &      minindx,maxindx,mindata,maxdata,                             &
     &      pminindx,pmaxindx,pmindata,pmaxdata,lon1,lon2,lat1,lat2,     &
     &      lev1,lev2,altlon1,altlon2,altlat1,altlat2,altlev1,altlev2,   &
     &      confidenceinterval,restrain,                                 &
     &      lag1,lag2,m1,m2,yr1,yr2,yr1a,yr2a,yr1b,yr2b,lsum,lsum2,lsel, &
     &      ndiff,ndiff2,indxuse,avex,avey,nfittime,minnum,intertype,    &
     &      mdiff,mdiff2,nens1,nens2,nyrwindow,lead1,lead2,irunvar,      &
     &      noisetype,debias,normalization,fitfunc,altavex,altavey,      &
     &      logscale,sqrtscale,lrank,ldetrend,lsubtract,lwrite,          &
     &      lincl,lks,dump,plot,lconting,fix2,lnormsd,lbootstrap,        &
     &      anom,composite,ensemble,squarescale,cubescale,               &
     &      twothirdscale,leads,                                         &
     &      logfield,sqrtfield,lfitnoise,lrandom,lrmse,                  &
     &      lchangesign,lmakeensfull,lweb,lstandardunits,lallobs,        &
     &      lstartstop,lnomissing,lensanom,lnooverlap,ncrossvalidate,    &
     &      day0,minfacsum,lmae,lincludelast,xyear,add_option,biasadd,   &
     &      biasmul,biasrt
        common /cc_getopts/ oper,strindx,plotfile,indexfiles,            &
     &      plotrunfile,corrmonths,namestring,bbfile
