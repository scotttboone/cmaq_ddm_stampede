MODEL = CCTM_ddm_stampede_Linux2_x86_64ifort

#   Compilers
 FC        = /opt/apps/intel13/mvapich2/1.9/bin/mpif90
 CC        = cc
 FPP       = /opt/apps/intel13/mvapich2/1.9/bin/mpif90

#lib_path  = $(lib_basedir)/${system}/${compiler}

 IOAPI_LIB = /home1/02554/stboone/apps/ioapi/Linux2_x86_64ifort
 PARIO_LIB = /work/02554/stboone/cmaq_4.7.1/lib/pario/Linux2_x86_64
 IOAPI_INC = /home1/02554/stboone/apps/ioapi/ioapi/fixed_src
 MPI_INC   = /opt/apps/intel13/mvapich2/1.9/include
 SE_INC    = /work/02554/stboone/cmaq_4.7.1/lib/stenex/Linux2_x86_64

#   Compiler flags
 add_flags = -traceback
 f_FLAGS   = ${add_flags} -extend_source 132 -vec-report0 -nodefines -cm -w95 -c -O2 -module /work/02554/stboone/cmaq_4.7.1/scripts/cctm/MOD_DIR -I $(SE_INC) -I $(MPI_INC) -I $(IOAPI_INC) -I $(IOAPI_LIB) -I $(PARIO_LIB) -I .
 F_FLAGS   = ${f_FLAGS}
f90_FLAGS  =
F90_FLAGS  = 
 C_FLAGS   = -v -g -DFLDMN -O0  -I /opt/apps/intel13/mvapich2/1.9/include
#
#   Linker and flags
 LINKER = ${FC}
 LINK_FLAGS =  -liomp5 -lpthread

#   Search PATH for source files
VPATH = .:ctm_yamo_ddm3d:par:init_yamo_ddm3d:gencoor_ddm3d:hyamo_ddm3d:vyamo_ddm3d:multiscale_ddm3d:acm2_inline:phot_table:ebi_cb05cl_ae5_ddm3d:aero5_ddm3d:aero_depv2:cloud_acm_ae5_ddm3d:pa:util:

#   Pre-compiler and flags
CPP = ${FC}
CPP_FLAGS = \
  -Dparallel \
  -Dcluster \
  -Dsens \
  -DINTERPB=PINTERPB \
  -DM3ERR=PM3ERR \
  -DM3EXIT=PM3EXIT \
  -DM3WARN=PM3WARN \
  -DSHUT3=PSHUT3 \
  -DWRITE3=PWRITE3 \
  -DSUBST_MODULES=SE_MODULES \
  -DSUBST_BARRIER=SE_BARRIER \
  -DSUBST_GLOBAL_MAX=SE_GLOBAL_MAX \
  -DSUBST_GLOBAL_MIN=SE_GLOBAL_MIN \
  -DSUBST_GLOBAL_MIN_DATA=SE_GLOBAL_MIN_DATA \
  -DSUBST_GLOBAL_TO_LOCAL_COORD=SE_GLOBAL_TO_LOCAL_COORD \
  -DSUBST_GLOBAL_SUM=SE_GLOBAL_SUM \
  -DSUBST_GLOBAL_LOGICAL=SE_GLOBAL_LOGICAL \
  -DSUBST_LOOP_INDEX=SE_LOOP_INDEX \
  -DSUBST_SUBGRID_INDEX=SE_SUBGRID_INDEX \
  -DSUBST_HI_LO_BND_PE=SE_HI_LO_BND_PE \
  -DSUBST_SUM_CHK=SE_SUM_CHK \
  -DSUBST_INIT_ARRAY=SE_INIT_ARRAY \
  -DSUBST_COMM=SE_COMM \
  -DSUBST_MY_REGION=SE_MY_REGION \
  -DSUBST_SLICE=SE_SLICE \
  -DSUBST_GATHER=SE_GATHER \
  -DSUBST_DATA_COPY=SE_DATA_COPY \
  -DSUBST_IN_SYN=SE_IN_SYN

#   Library paths
#LIBRARIES = \
#-L$(lib_path)/pario -lpario \
#-L$(lib_path)/ioapi_3.1/$(LIOAPI)_50000 -lioapi \
#-L$(lib_path)/netcdf/lib -lnetcdf \
#-L$(lib_path)/se_snl -lse_snl \
#-L$(lib_path)/mpich/lib -lmpich $(extra_lib)

 SE_SNL    = -L/work/02554/stboone/cmaq_4.7.1/lib/stenex/Linux2_x86_64 -lse_snl
 PARIO     = -L/work/02554/stboone/cmaq_4.7.1/lib/pario/Linux2_x86_64 -lpario
#MPICH     = -L/usr/intel/mvapich-1.2.0/lib -lmpich
 IOAPI     = -L/home1/02554/stboone/apps/ioapi/Linux2_x86_64ifort -lioapi
 NETCDF    = -L/home1/02554/stboone/apps/netcdf/lib -lnetcdf
 LIBRARIES = $(SE_SNL) $(PARIO) $(MPICH) $(IOAPI) $(NETCDF)

# Includes paths
 BASE_INC = /work/02554/stboone/cmaq_4.7.1/scripts/cctm/BLD
 MECH_INC = $(BASE_INC)
 TRAC_INC = $(BASE_INC)

#   Includes
INCLUDES = \
  -DSUBST_GRID_ID= \
  -DSUBST_PE_COMM=\"$(BASE_INC)/PE_COMM.EXT\" \
  -DSUBST_CONST=\"$(BASE_INC)/CONST.EXT\" \
  -DSUBST_FILES_ID=\"$(BASE_INC)/FILES_CTM.EXT\" \
  -DSUBST_EMPR_VD=\"$(BASE_INC)/EMISPRM.vdif.EXT\" \
  -DSUBST_EMPR_CH=\"$(BASE_INC)/EMISPRM.chem.EXT\" \
  -DSUBST_IOPARMS=\"$(IOAPI_INC)/PARMS3.EXT\" \
  -DSUBST_IOFDESC=\"$(IOAPI_INC)/FDESC3.EXT\" \
  -DSUBST_IODECL=\"$(IOAPI_INC)/IODECL3.EXT\" \
  -DSUBST_RXCMMN=\"$(BASE_INC)/RXCM.EXT\" \
  -DSUBST_RXDATA=\"$(BASE_INC)/RXDT.EXT\" \
  -DSUBST_GC_SPC=\"$(MECH_INC)/GC_SPC.EXT\" \
  -DSUBST_GC_EMIS=\"$(MECH_INC)/GC_EMIS.EXT\" \
  -DSUBST_GC_ICBC=\"$(MECH_INC)/GC_ICBC.EXT\" \
  -DSUBST_GC_DIFF=\"$(MECH_INC)/GC_DIFF.EXT\" \
  -DSUBST_GC_DDEP=\"$(MECH_INC)/GC_DDEP.EXT\" \
  -DSUBST_GC_DEPV=\"$(MECH_INC)/GC_DEPV.EXT\" \
  -DSUBST_GC_ADV=\"$(MECH_INC)/GC_ADV.EXT\" \
  -DSUBST_GC_CONC=\"$(MECH_INC)/GC_CONC.EXT\" \
  -DSUBST_GC_G2AE=\"$(MECH_INC)/GC_G2AE.EXT\" \
  -DSUBST_GC_G2AQ=\"$(MECH_INC)/GC_G2AQ.EXT\" \
  -DSUBST_GC_SCAV=\"$(MECH_INC)/GC_SCAV.EXT\" \
  -DSUBST_GC_WDEP=\"$(MECH_INC)/GC_WDEP.EXT\" \
  -DSUBST_AE_SPC=\"$(MECH_INC)/AE_SPC.EXT\" \
  -DSUBST_AE_EMIS=\"$(MECH_INC)/AE_EMIS.EXT\" \
  -DSUBST_AE_ICBC=\"$(MECH_INC)/AE_ICBC.EXT\" \
  -DSUBST_AE_DIFF=\"$(MECH_INC)/AE_DIFF.EXT\" \
  -DSUBST_AE_DDEP=\"$(MECH_INC)/AE_DDEP.EXT\" \
  -DSUBST_AE_DEPV=\"$(MECH_INC)/AE_DEPV.EXT\" \
  -DSUBST_AE_ADV=\"$(MECH_INC)/AE_ADV.EXT\" \
  -DSUBST_AE_CONC=\"$(MECH_INC)/AE_CONC.EXT\" \
  -DSUBST_AE_A2AQ=\"$(MECH_INC)/AE_A2AQ.EXT\" \
  -DSUBST_AE_SCAV=\"$(MECH_INC)/AE_SCAV.EXT\" \
  -DSUBST_AE_WDEP=\"$(MECH_INC)/AE_WDEP.EXT\" \
  -DSUBST_NR_SPC=\"$(MECH_INC)/NR_SPC.EXT\" \
  -DSUBST_NR_EMIS=\"$(MECH_INC)/NR_EMIS.EXT\" \
  -DSUBST_NR_ICBC=\"$(MECH_INC)/NR_ICBC.EXT\" \
  -DSUBST_NR_DIFF=\"$(MECH_INC)/NR_DIFF.EXT\" \
  -DSUBST_NR_DDEP=\"$(MECH_INC)/NR_DDEP.EXT\" \
  -DSUBST_NR_DEPV=\"$(MECH_INC)/NR_DEPV.EXT\" \
  -DSUBST_NR_ADV=\"$(MECH_INC)/NR_ADV.EXT\" \
  -DSUBST_NR_N2AE=\"$(MECH_INC)/NR_N2AE.EXT\" \
  -DSUBST_NR_CONC=\"$(MECH_INC)/NR_CONC.EXT\" \
  -DSUBST_NR_N2AQ=\"$(MECH_INC)/NR_N2AQ.EXT\" \
  -DSUBST_NR_SCAV=\"$(MECH_INC)/NR_SCAV.EXT\" \
  -DSUBST_NR_WDEP=\"$(MECH_INC)/NR_WDEP.EXT\" \
  -DSUBST_TR_SPC=\"$(TRAC_INC)/TR_SPC.EXT\" \
  -DSUBST_TR_EMIS=\"$(TRAC_INC)/TR_EMIS.EXT\" \
  -DSUBST_TR_ICBC=\"$(TRAC_INC)/TR_ICBC.EXT\" \
  -DSUBST_TR_DIFF=\"$(TRAC_INC)/TR_DIFF.EXT\" \
  -DSUBST_TR_DDEP=\"$(TRAC_INC)/TR_DDEP.EXT\" \
  -DSUBST_TR_DEPV=\"$(TRAC_INC)/TR_DEPV.EXT\" \
  -DSUBST_TR_ADV=\"$(TRAC_INC)/TR_ADV.EXT\" \
  -DSUBST_TR_T2AQ=\"$(TRAC_INC)/TR_T2AQ.EXT\" \
  -DSUBST_TR_SCAV=\"$(TRAC_INC)/TR_SCAV.EXT\" \
  -DSUBST_TR_WDEP=\"$(TRAC_INC)/TR_WDEP.EXT\" \
  -DSUBST_PACTL_ID=\"/work/02554/stboone/cmaq_4.7.1/models/pa_noop/PA_CTL.EXT\" \
  -DSUBST_PACMN_ID=\"/work/02554/stboone/cmaq_4.7.1/models/pa_noop/PA_CMN.EXT\" \
  -DSUBST_PADAT_ID=\"/work/02554/stboone/cmaq_4.7.1/models/pa_noop/PA_DAT.EXT\" \
  -DSUBST_MPICH=\"/opt/apps/intel13/mvapich2/1.9/include/mpif.h\"

GLOBAL_MODULES = \
  CGRID_SPCS.o \
  HGRD_DEFN.o \
  VGRD_DEFN.o \
  STD_CONC.o \
  GRID_CONF.o \
  WVEL_DEFN.o \
  PMEM_DEFN.o \
  UDTYPES.o \
  STK_EMIS.o \
  STK_PRMS.o \
  PTMAP.o \
  PTBILIN.o \
  PTMET.o \
  PT3D_DEFN.o \
  AERO_EMIS.o \
  PAGRD_DEFN.o

CTM_YAMO_DDM3D = \
  AVG_CONC.o \
  PCGRID_DEFN.o \
  S_PCGRID_DEFN.o \
  advstep.o \
  driver.o \
  hcontvel.o \
  hveloc.o \
  sciproc.o \
  wr_aconc.o \
  wr_asens.o \
  wr_cgrid.o \
  wr_conc.o \
  wr_sengrid.o

PAR_NODISTR = \
  distr_env.o \
  par_init.o \
  par_term.o 

INIT_YAMO_DDM3D = \
  cktime.o \
  diffmsg.o \
  flcheck.o \
  grdcheck.o \
  initscen.o \
  load_cgrid.o \
  load_sengrid.o \
  opconc.o \
  s_lu.o \
  setup_logdev.o \
  sinput.o

GENCOOR_DDM3D = \
  conv_sengrid.o \
  couple.o \
  ppmv_msmr.o

HYAMO_DDM3D = \
  advbc_map.o \
  hadvyppm.o \
  hppm.o \
  rdbcon.o \
  s_zfdbc.o \
  x_ppm.o \
  x_yamo.o \
  y_ppm.o \
  y_yamo.o \
  zfdbc.o

VYAMO_DDM3D = \
  vppm.o \
  zadvyppm.o

MULTISCALE_DDM3D = \
  deform.o \
  hcdiff3d.o \
  hdiff.o \
  rho_j.o

ACM2_INLINE = \
  BIOG_EMIS.o \
  BEIS_DEFN.o \
  DEPVVARS.o \
  DEPV_DEFN.o \
  EMIS_DEFN.o \
  VDIFF_MAP.o \
  S_EMIS_DEFN.o \
  beis3.o \
  cgrid_depv.o \
  checkmem.o \
  chkgrid.o \
  conv_cgrid.o \
  czangle.o \
  delta_zs.o \
  eddyx.o \
  fire_plmris.o \
  gas_depv_map.o \
  getfline.o \
  getparb.o \
  hrno.o \
  m3dry.o \
  matrix.o \
  smatrix.o \
  opddep.o \
  opdepv_diag.o \
  opemis.o \
  s_opemis.o \
  openlayout.o \
  opmet.o \
  oppt3d_diag.o \
  opssemis.o \
  parsline.o \
  plmris.o \
  plsprd.o \
  preplm.o \
  rddepv.o \
  rdmet.o \
  ssemis.o \
  tmpbeis.o \
  tri.o \
  stri.o \
  ungridb2.o \
  vdiffacm2.o \
  wrdaymsg.o \
  write3_distr.o

PHOT_TABLE = \
  phot.o

EBI_CB05CL_AE5_DDM3D = \
  hrdata_mod.o \
  hrcalcks.o \
  hrdriver.o \
  hrg1.o \
  hrg2.o \
  hrg3.o \
  hrg4.o \
  hrinit.o \
  hrprodloss.o \
  hrrates.o \
  hrsolver.o \
  jac_cb05cl_ae5_aq.o

AERO5_DDM3D = \
  AERO_INFO.o \
  SOA_NEWT.o \
  aero_driver.o \
  aero_sens.o \
  aero_subs.o \
  hetchem.o \
  isocom.o \
  isofwd.o \
  isorev.o \
  opdiam.o \
  opvis.o \
  orgaer5.o

AERO_DEPV2 = \
  aero_depv.o

CLOUD_ACM_AE5_DDM3D = \
  acmcld.o \
  aq_map.o \
  aqchem.o \
  cldproc_acm.o \
  convcld_acm.o \
  getalpha.o \
  hlconst.o \
  indexn.o \
  rescld.o \
  scavwdep.o \
  set_aeconcmin.o

PA = \
  pa_init.o \
  pa_irr.o \
  pa_irr_ctl.o \
  pa_mkhdr.o \
  pa_output.o \
  pa_update.o

UTIL = \
  cksummer.o \
  findex.o \
  get_envlist.o \
  lstepf.o \
  subhdomain.o \
  subhfile.o

OBJS = \
  $(GLOBAL_MODULES) \
  $(CTM_YAMO_DDM3D) \
  $(PAR_NODISTR) \
  $(INIT_YAMO_DDM3D) \
  $(GENCOOR_DDM3D) \
  $(HYAMO_DDM3D) \
  $(VYAMO_DDM3D) \
  $(MULTISCALE_DDM3D) \
  $(ACM2_INLINE) \
  $(PHOT_TABLE) \
  $(EBI_CB05CL_AE5_DDM3D) \
  $(AERO5_DDM3D) \
  $(AERO_DEPV2) \
  $(CLOUD_ACM_AE5_DDM3D) \
  $(PA) \
  $(UTIL)

.SUFFIXES: .F .f .c .F90 .f90

$(MODEL): $(OBJS)
	$(LINKER) $(LINKER_FLAGS) $(OBJS) $(LIBRARIES) -o $@

.F.o:
	$(FC) -c $(F_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f.o:
	$(FC) -c $(f_FLAGS) $<

.F90.o:
	$(FC) -c $(F90_FLAGS) $(CPP_FLAGS) $(INCLUDES) $<

.f90.o:
	$(FC) -c $(f90_FLAGS) $<

.c.o:
	$(CC) -c $(C_FLAGS) $<

clean:
	rm -f $(OBJS) $(MODEL) *.mod

