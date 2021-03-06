#include "cppdefs.h"

      subroutine check_srcs

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if the new files     !!!!!!
!!!!!! were introduced into or deleted from the list of source !!!!!!
!!!!!! codes SRCS in the Makefile.                             !!!!!!

      implicit none
      integer i
#include "strings.h"
      do i=1,max_opt_size
        srcs(i:i)=' '
      enddo
      srcs(   1:   2)='90'
      srcs(   3:   5)='  \'
      srcs(   6:  19)=' module_nh.F90'
      srcs(  20:  34)=' module_nbq.F90'
      srcs(  35:  58)=' module_tracetxt_out.F90'
      srcs(  59:  70)=' nump_nh.F90'
      srcs(  71:  91)=' mat_cont_init_nh.F90'
      srcs(  92: 107)=' mat_cont_nh.F90'
      srcs( 108: 127)=' mat_mom_init_nh.F90'
      srcs( 128: 142)=' mat_mom_nh.F90'
      srcs( 143: 157)=' initial_nh.F90'
      srcs( 158: 176)=' z_thickness_nh.F90'
      srcs( 177: 200)=' module_tracetxt_out.F90'
      srcs( 201: 211)=' ru_nbq.F90'
      srcs( 212: 227)=' viscous_nbq.F90'
      srcs( 228: 242)=' step3d_nbq.F90'
      srcs( 243: 258)=' density_nbq.F90'
      srcs( 259: 276)=' spkitlocal_nh.F90'
      srcs( 277: 293)=' module_param.F90'
      srcs( 294: 300)=' main.F'
      srcs( 301: 307)=' step.F'
      srcs( 308: 318)=' read_inp.F'
      srcs( 319: 332)=' timers_roms.F'
      srcs( 333: 347)=' init_scalars.F'
      srcs( 348: 361)=' init_arrays.F'
      srcs( 362: 375)=' set_weights.F'
      srcs( 376: 388)=' set_scoord.F'
      srcs( 389: 399)=' ana_grid.F'
      srcs( 400: 413)=' setup_grid1.F'
      srcs( 414: 427)=' setup_grid2.F'
      srcs( 428: 441)=' set_nudgcof.F'
      srcs( 442: 455)=' ana_initial.F'
      srcs( 456: 468)=' analytical.F'
      srcs( 469: 477)=' zonavg.F'
      srcs( 478: 490)=' setup_kwds.F'
      srcs( 491: 503)=' check_kwds.F'
      srcs( 504: 516)=' check_srcs.F'
      srcs( 517: 534)=' check_switches1.F'
      srcs( 535: 552)=' check_switches2.F'
      srcs( 553: 563)=' step3d_w.F'
      srcs( 564: 576)=' rhs3d_w_nh.F'
      srcs( 577: 584)=' w3dbc.F'
      srcs( 585: 603)=' cpl_prism_define.F'
      srcs( 604: 619)=' cpl_prism_put.F'
      srcs( 620: 636)=' cpl_prism_init.F'
      srcs( 637: 652)=' cpl_prism_get.F'
      srcs( 653: 671)=' cpl_prism_getvar.F'
      srcs( 672: 695)=' cpl_prism_getlistvals.F'
      srcs( 696: 704)=' step2d.F'
      srcs( 705: 712)=' u2dbc.F'
      srcs( 713: 720)=' v2dbc.F'
      srcs( 721: 729)=' zetabc.F'
      srcs( 730: 743)=' obc_volcons.F'
      srcs( 744: 756)=' pre_step3d.F'
      srcs( 757: 767)=' step3d_t.F'
      srcs( 768: 780)=' step3d_uv1.F'
      srcs( 781: 793)=' step3d_uv2.F'
      srcs( 794: 802)=' prsgrd.F'
      srcs( 803: 810)=' rhs3d.F'
      srcs( 811: 822)=' set_depth.F'
      srcs( 823: 830)=' omega.F'
      srcs( 831: 840)=' uv3dmix.F'
      srcs( 841: 854)=' uv3dmix_spg.F'
      srcs( 855: 863)=' t3dmix.F'
      srcs( 864: 876)=' t3dmix_spg.F'
      srcs( 877: 888)=' hmix_coef.F'
      srcs( 889: 897)=' wetdry.F'
      srcs( 898: 905)=' u3dbc.F'
      srcs( 906: 913)=' v3dbc.F'
      srcs( 914: 921)=' t3dbc.F'
      srcs( 922: 931)=' rho_eos.F'
      srcs( 932: 942)=' ab_ratio.F'
      srcs( 943: 953)=' alfabeta.F'
      srcs( 954: 964)=' ana_vmix.F'
      srcs( 965: 974)=' bvf_mix.F'
      srcs( 975: 985)=' lmd_vmix.F'
      srcs( 986: 998)=' gls_mixing.F'
      srcs( 999:1009)=' lmd_skpp.F'
      srcs(1010:1020)=' lmd_bkpp.F'
      srcs(1021:1033)=' lmd_swfrac.F'
      srcs(1034:1046)=' lmd_wscale.F'
      srcs(1047:1060)=' gls_prestep.F'
      srcs(1061:1067)=' diag.F'
      srcs(1068:1076)=' wvlcty.F'
      srcs(1077:1088)=' checkdims.F'
      srcs(1089:1105)=' grid_stiffness.F'
      srcs(1106:1116)=' bio_diag.F'
      srcs(1117:1127)=' get_date.F'
      srcs(1128:1136)=' lenstr.F'
      srcs(1137:1147)=' closecdf.F'
      srcs(1148:1163)=' ana_initracer.F'
      srcs(1164:1177)=' insert_node.F'
      srcs(1178:1188)=' nf_fread.F'
      srcs(1189:1199)=' get_grid.F'
      srcs(1200:1213)=' get_initial.F'
      srcs(1214:1227)=' def_grid_2d.F'
      srcs(1228:1241)=' def_grid_3d.F'
      srcs(1242:1251)=' def_his.F'
      srcs(1252:1261)=' def_rst.F'
      srcs(1262:1271)=' set_avg.F'
      srcs(1272:1282)=' wrt_grid.F'
      srcs(1283:1292)=' wrt_his.F'
      srcs(1293:1302)=' wrt_rst.F'
      srcs(1303:1312)=' wrt_avg.F'
      srcs(1313:1321)=' output.F'
      srcs(1322:1339)=' put_global_atts.F'
      srcs(1340:1352)=' nf_fread_x.F'
      srcs(1353:1365)=' nf_fread_y.F'
      srcs(1366:1377)=' def_diags.F'
      srcs(1378:1389)=' wrt_diags.F'
      srcs(1390:1405)=' wrt_diags_avg.F'
      srcs(1406:1421)=' set_diags_avg.F'
      srcs(1422:1434)=' def_diagsM.F'
      srcs(1435:1447)=' wrt_diagsM.F'
      srcs(1448:1464)=' wrt_diagsM_avg.F'
      srcs(1465:1481)=' set_diagsM_avg.F'
      srcs(1482:1497)=' def_bio_diags.F'
      srcs(1498:1513)=' wrt_bio_diags.F'
      srcs(1514:1533)=' wrt_bio_diags_avg.F'
      srcs(1534:1553)=' set_bio_diags_avg.F'
      srcs(1554:1563)=' get_vbc.F'
      srcs(1564:1575)=' set_cycle.F'
      srcs(1576:1587)=' get_wwave.F'
      srcs(1588:1600)=' get_tclima.F'
      srcs(1601:1613)=' get_uclima.F'
      srcs(1614:1623)=' get_ssh.F'
      srcs(1624:1633)=' get_sss.F'
      srcs(1634:1646)=' get_smflux.F'
      srcs(1647:1659)=' get_stflux.F'
      srcs(1660:1672)=' get_srflux.F'
      srcs(1673:1682)=' get_sst.F'
      srcs(1683:1694)=' get_tides.F'
      srcs(1695:1706)=' clm_tides.F'
      srcs(1707:1717)=' get_bulk.F'
      srcs(1718:1729)=' bulk_flux.F'
      srcs(1730:1739)=' get_bry.F'
      srcs(1740:1753)=' get_bry_bio.F'
      srcs(1754:1767)=' nf_read_bry.F'
      srcs(1768:1777)=' sstskin.F'
      srcs(1778:1791)=' get_psource.F'
      srcs(1792:1808)=' get_psource_ts.F'
      srcs(1809:1827)=' nf_add_attribute.F'
      srcs(1828:1839)=' fillvalue.F'
      srcs(1840:1849)=' mrl_wci.F'
      srcs(1850:1861)=' wkb_wwave.F'
      srcs(1862:1869)=' wkbbc.F'
      srcs(1870:1883)=' get_bry_wkb.F'
      srcs(1884:1901)=' online_bulk_var.F'
      srcs(1902:1919)=' online_get_bulk.F'
      srcs(1920:1935)=' online_interp.F'
      srcs(1936:1961)=' online_interpolate_bulk.F'
      srcs(1962:1979)=' online_set_bulk.F'
      srcs(1980:1991)=' MPI_Setup.F'
      srcs(1992:2004)=' MessPass2D.F'
      srcs(2005:2017)=' MessPass3D.F'
      srcs(2018:2028)=' exchange.F'
      srcs(2029:2038)=' biology.F'
      srcs(2039:2047)=' o2sato.F'
      srcs(2048:2058)=' sediment.F'
      srcs(2059:2064)=' bbl.F'
      srcs(2065:2078)=' init_floats.F'
      srcs(2079:2091)=' wrt_floats.F'
      srcs(2092:2105)=' step_floats.F'
      srcs(2106:2118)=' rhs_floats.F'
      srcs(2119:2131)=' interp_rho.F'
      srcs(2132:2144)=' def_floats.F'
      srcs(2145:2165)=' init_arrays_floats.F'
      srcs(2166:2179)=' random_walk.F'
      srcs(2180:2200)=' get_initial_floats.F'
      srcs(2201:2211)=' init_sta.F'
      srcs(2212:2221)=' wrt_sta.F'
      srcs(2222:2232)=' step_sta.F'
      srcs(2233:2245)=' interp_sta.F'
      srcs(2246:2255)=' def_sta.F'
      srcs(2256:2273)=' init_arrays_sta.F'
      srcs(2274:2280)=' zoom.F'
      srcs(2281:2291)=' update2D.F'
      srcs(2292:2310)=' set_nudgcof_fine.F'
      srcs(2311:2322)=' zoombc_2D.F'
      srcs(2323:2334)=' zoombc_3D.F'
      srcs(2335:2347)=' uv3dpremix.F'
      srcs(2348:2359)=' t3dpremix.F'
      srcs(2360:2370)=' update3D.F'
      srcs(2371:2381)=' ropigrid.F'
      srcs(2382:2392)=' ropistep.F'
      srcs(2393:2403)=' pirostep.F'
      srcs(2404:2412)=' trcini.F'
      srcs(2413:2421)=' trclsm.F'
      srcs(2422:2434)=' pisces_ini.F'
      srcs(2435:2443)=' p4zday.F'
      srcs(2444:2453)=' p4znano.F'
      srcs(2454:2462)=' p4zlys.F'
      srcs(2463:2472)=' p4zdiat.F'
      srcs(2473:2481)=' p4zopt.F'
      srcs(2482:2491)=' p4zsink.F'
      srcs(2492:2500)=' p4zflx.F'
      srcs(2501:2509)=' p4zprg.F'
      srcs(2510:2519)=' p4zslow.F'
      srcs(2520:2528)=' p4zint.F'
      srcs(2529:2538)=' p4zprod.F'
      srcs(2539:2549)=' p4zmicro.F'
      srcs(2550:2558)=' p4zlim.F'
      srcs(2559:2567)=' p4zrem.F'
      srcs(2568:2576)=' p4zche.F'
      srcs(2577:2585)=' p4zbio.F'
      srcs(2586:2595)=' p4zmeso.F'
      srcs(2596:2606)=' p4zsink2.F'
      srcs(2607:2615)=' p4zsed.F'
      srcs(2616:2628)=' autotiling.F'
      srcs(2629:2636)=' debug.F'
      srcs(2637:2650)=' Agrif2Model.F'
      return
      end
