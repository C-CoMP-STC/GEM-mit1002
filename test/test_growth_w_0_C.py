import unittest

import cobra

class TestGrowthWithNoCarbon(unittest.TestCase):
    @unittest.expectedFailure # Fails because the model is so small that it is missing some of these metabolites
    # Hopefully will work fine once I update the model to actually be MIT1002
    def test_growth_w_0_C(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model('../model.xml')

        # Set the media so that there are no carbon sources
        medium = {'EX_12dgr160_e': 0.0,
                'EX_12dgr180_e': 0.0,
                'EX_1ag160_e': 0.0,
                'EX_1ag180_e': 0.0,
                'EX_1ag181d9_e': 0.0,
                'EX_1ag182d9d12_e': 0.0,
                'EX_25dkglcn_e': 0.0,
                'EX_2ameph_e': 0.0, # This wan't in BiGG?
                'EX_2ddglcn_e': 0.0,
                'EX_2m35mdntha_e': 0.0,
                'EX_2pglyc_e': 0.0,
                'EX_34dhbz_e': 0.0,
                'EX_35dnta_e': 0.0,
                'EX_3h4atb_e': 0.0,
                'EX_3mb_e': 0.0,
                'EX_3oxoadp_e': 0.0,
                'EX_4abut_e': 0.0,
                'EX_4hbald_e': 0.0,
                'EX_4hbz_e': 0.0,
                'EX_5aptn_e': 0.0,
                'EX_5drib_e': 0.0,
                'EX_6atha_e': 0.0,
                'EX_6hnac_e': 0.0,
                'EX_LalaLglu_e': 0.0,
                'EX_Lcyst_e': 0.0,
                'EX_R_3hphpa_e': 0.0,
                'EX_R_3hppta_e': 0.0,
                'EX_ac_e': 0.0, # Acetate
                'EX_acac_e': 0.0,
                'EX_acald_e': 0.0,
                'EX_acgam_e': 0.0,
                'EX_acnam_e': 0.0,
                'EX_acon_C_e': 0.0,
                'EX_actn__R_e': 0.0,
                'EX_ad_e': 0.0,
                'EX_ade_e': 0.0,
                'EX_adn_e': 0.0,
                'EX_akg_e': 0.0,
                'EX_ala_B_e': 0.0,
                'EX_ala__L_e': 0.0,
                'EX_alaala_e': 0.0,
                'EX_alahis_e': 0.0,
                'EX_alaleu_e': 0.0,
                'EX_alathr_e': 0.0,
                'EX_alatrp_e': 0.0,
                'EX_algac_MG_14_e': 0.0,
                'EX_algac_MG_23_e': 0.0,
                'EX_algac_MG_32_e': 0.0,
                'EX_algac_MG_41_e': 0.0,
                'EX_algac__M_e': 0.0,
                'EX_alltn_e': 0.0,
                'EX_anhgm_e': 0.0,
                'EX_arab__L_e': 0.0,
                'EX_arbt6p_e': 0.0,
                'EX_arbt_e': 0.0,
                'EX_arg__L_e': 0.0,
                'EX_asn__L_e': 0.0,
                'EX_aso3_e': 1000.0,
                'EX_aso4_e': 1000.0,
                'EX_asp__L_e': 0.0,
                'EX_balaala_e': 0.0,
                'EX_balabala_e': 0.0,
                'EX_balagly_e': 0.0,
                'EX_balaleu_e': 0.0,
                'EX_balamd_e': 0.0,
                'EX_but_e': 0.0,
                'EX_butso3_e': 0.0,
                'EX_bz_e': 0.0,
                'EX_ca2_e': 1000.0,
                'EX_carn_e': 0.0,
                'EX_cell4_e': 0.0,
                'EX_cellb_e': 0.0,
                'EX_cgly_e': 0.0,
                'EX_chol_e': 0.0,
                'EX_chols_e': 0.0,
                'EX_chsterol_e': 0.0,
                'EX_cinnm_e': 0.0,
                'EX_cit_e': 0.0,
                'EX_cl_e': 1000.0,
                'EX_cmcbtt_e': 0.0,
                'EX_co2_e': 0.0, # CO2
                'EX_co_e': 0.0,
                'EX_cobalt2_e': 1000.0,
                'EX_confrl_e': 0.0,
                'EX_crn_e': 0.0,
                'EX_csn_e': 0.0,
                'EX_cu2_e': 1000.0,
                'EX_cyan_e': 0.0, # Hydrogen cyanide (HCN)
                'EX_cynt_e': 0.0, # Cyanate (CNA)
                'EX_cys__L_e': 0.0,
                'EX_d23hb_e': 0.0,
                'EX_dad_2_e': 0.0,
                'EX_dag181d9_e': 0.0,
                'EX_dag182d9d12_e': 0.0,
                'EX_dca_e': 0.0,
                'EX_ddca_e': 0.0,
                'EX_dmanur_e': 0.0,
                'EX_dmgly_e': 0.0,
                'EX_dmso2_e': 0.0,
                'EX_drib_e': 0.0,
                'EX_dxyl_e': 0.0,
                'EX_ecto__L_e': 0.0,
                'EX_enter_e': 0.0,
                'EX_etha_e': 0.0,
                'EX_ethso3_e': 0.0,
                'EX_etoh_e': 0.0,
                'EX_fald_e': 0.0,
                'EX_fcmcbtt_e': 0.0,
                'EX_fe2_e': 1000.0,
                'EX_fe3_e': 1000.0,
                'EX_fe3dcit_e': 0.0,
                'EX_fe3dhbzs3_e': 0.0,
                'EX_fe3pyovd_kt_e': 0.0,
                'EX_feenter_e': 0.0,
                'EX_feoxam_e': 0.0,
                'EX_feoxam_un_e': 0.0,
                'EX_fer_e': 0.0,
                'EX_for_e': 0.0,
                'EX_fru_e': 0.0,
                'EX_fuc_e': 0.0,
                'EX_fum_e': 0.0,
                'EX_g3pc_e': 0.0,
                'EX_g3pg_e': 0.0,
                'EX_g3pi_e': 0.0,
                'EX_g3ps_e': 0.0,
                'EX_gal_bD_e': 0.0,
                'EX_gal_e': 0.0,
                'EX_galct__D_e': 0.0,
                'EX_galctn__D_e': 0.0,
                'EX_galctn__L_e': 0.0,
                'EX_galctr__D_e': 0.0,
                'EX_galt_e': 0.0,
                'EX_galur_e': 0.0,
                'EX_gam_e': 0.0,
                'EX_glc__D_e': 0.0, # Glucose- now set to 0
                'EX_glc__aD_e': 0.0, # Alpha glucose
                'EX_glcn_e': 0.0,
                'EX_glcr_e': 0.0,
                'EX_glcur_e': 0.0,
                'EX_gln__L_e': 0.0,
                'EX_glu__L_e': 0.0,
                'EX_glutar_e': 0.0,
                'EX_gly_e': 0.0,
                'EX_glyald_e': 0.0,
                'EX_glyb_e': 0.0,
                'EX_glyc2p_e': 0.0,
                'EX_glyc3p_e': 0.0,
                'EX_glyc__R_e': 0.0,
                'EX_glyc_e': 0.0,
                'EX_glyclt_e': 0.0,
                'EX_glycol_e': 0.0,
                'EX_glygln_e': 0.0,
                'EX_glyglu_e': 0.0,
                'EX_glygly_e': 0.0,
                'EX_glymet_e': 0.0,
                'EX_glyphe_e': 0.0,
                'EX_glyser_e': 0.0,
                'EX_gsn_e': 0.0,
                'EX_gthrd_e': 0.0,
                'EX_gua_e': 0.0,
                'EX_h2_e': 1000.0,
                'EX_h2o2_e': 1000.0,
                'EX_h2o_e': 1000.0,
                'EX_h2s_e': 1000.0,
                'EX_h_e': 1000.0,
                'EX_hco3_e': 0.0, # Bicarbonate
                'EX_hdca_e': 0.0,
                'EX_hdcea_e': 0.0,
                'EX_hexs_e': 0.0,
                'EX_hia_e': 0.0, # Pseudoreaction
                'EX_his__L_e': 0.0,
                'EX_hisgly_e': 0.0,
                'EX_hishis_e': 0.0,
                'EX_hom__L_e': 0.0,
                'EX_hpta_e': 0.0,
                'EX_hxa_e': 0.0,
                'EX_hxan_e': 0.0,
                'EX_icit_e': 0.0,
                'EX_ile__L_e': 0.0,
                'EX_ind3ac_e': 0.0,
                'EX_indole_e': 0.0,
                'EX_inost_e': 0.0,
                'EX_ins_e': 0.0,
                'EX_isetac_e': 0.0,
                'EX_istfrnA_e': 0.0,
                'EX_istfrnB_e': 0.0,
                'EX_k_e': 1000.0,
                'EX_lac__D_e': 0.0,
                'EX_lac__L_e': 0.0,
                'EX_lcts_e': 0.0,
                'EX_leu__L_e': 0.0,
                'EX_leuleu_e': 0.0,
                'EX_lnlc_e': 0.0,
                'EX_lys__D_e': 0.0,
                'EX_lys__L_e': 0.0,
                'EX_m_xyl_e': 0.0,
                'EX_mal__D_e': 0.0,
                'EX_mal__L_e': 0.0,
                'EX_malt_e': 0.0,
                'EX_malthp_e': 0.0, # Pseudoreaction
                'EX_malthx_e': 0.0,
                'EX_malttr_e': 0.0,
                'EX_maltttr_e': 0.0,
                'EX_man_e': 0.0,
                'EX_manur_e': 0.0,
                'EX_melib_e': 0.0,
                'EX_meoh_e': 0.0,
                'EX_met__D_e': 0.0,
                'EX_met__L_e': 0.0,
                'EX_metox__R_e': 0.0,
                'EX_metox_e': 0.0,
                'EX_mg2_e': 1000.0,
                'EX_mn2_e': 1000.0,
                'EX_mnl_e': 0.0,
                'EX_mobd_e': 1000.0,
                'EX_mso3_e': 0.0,
                'EX_nac_e': 0.0,
                'EX_nh4_e': 1000.0,
                'EX_nmn_e': 0.0,
                'EX_no2_e': 1000.0,
                'EX_no3_e': 1000.0,
                'EX_no_e': 1000.0,
                'EX_nona_e': 0.0,
                'EX_o2_e': 1000.0,
                'EX_oaa_e': 0.0,
                'EX_ocdca_e': 0.0,
                'EX_ocdcea_e': 0.0,
                'EX_octa_e': 0.0,
                'EX_orn__D_e': 0.0,
                'EX_orn_e': 0.0,
                'EX_orot_e': 0.0,
                'EX_oxa_e': 0.0,
                'EX_p_xyl_e': 0.0, # pseudoreaction
                'EX_pac_e': 0.0,
                'EX_pacald_e': 0.0,
                'EX_pea_e': 0.0,
                'EX_peamn_e': 0.0,
                'EX_pentso3_e': 0.0,
                'EX_phe__L_e': 0.0,
                'EX_phedca_e': 0.0,
                'EX_phehpa_e': 0.0,
                'EX_phehxa_e': 0.0,
                'EX_pheme_e': 0.0,
                'EX_phenona_e': 0.0,
                'EX_pheocta_e': 0.0,
                'EX_phept_e': 0.0,
                'EX_pi_e': 1000.0,
                'EX_pnto__R_e': 0.0,
                'EX_ppa_e': 0.0,
                'EX_ppi_e': 1000.0, # Diphosphate
                'EX_prealg_MG_14_e': 0.0,
                'EX_prealg_MG_23_e': 0.0,
                'EX_prealg_MG_32_e': 0.0,
                'EX_prealg_MG_41_e': 0.0,
                'EX_prealginate_G_e': 0.0,
                'EX_pro__L_e': 0.0,
                'EX_progly_e': 0.0,
                'EX_pta_e': 0.0,
                'EX_ptrc_e': 0.0,
                'EX_ptsla_e': 0.0,
                'EX_pydxn_e': 0.0,
                'EX_pyovd_kt_e': 0.0,
                'EX_pyr_e': 0.0, # Pyruvate
                'EX_quin_e': 0.0,
                'EX_raffin_e': 0.0,
                'EX_rib__D_e': 0.0,
                'EX_rnam_e': 0.0,
                'EX_salchs2_e': 0.0,
                'EX_salchs4fe_e': 0.0,
                'EX_salchsx_e': 0.0,
                'EX_sbt__D_e': 0.0,
                'EX_sel_e': 1000.0,
                'EX_ser__L_e': 0.0,
                'EX_sheme_e': 0.0,
                'EX_skm_e': 0.0,
                'EX_slnt_e': 1000.0,
                'EX_so3_e': 1000.0,
                'EX_so4_e': 1000.0,
                'EX_spmd_e': 0.0,
                'EX_stfrnA_e': 0.0,
                'EX_stfrnB_e': 0.0,
                'EX_succ_e': 0.0, # Succinate
                'EX_sucr_e': 0.0, # Sucrose
                'EX_sulfac_e': 0.0,
                'EX_tag160_e': 0.0,
                'EX_tag180_e': 0.0,
                'EX_tag181d9_e': 0.0,
                'EX_tag182d9d12_e': 0.0,
                'EX_tartr__L_e': 0.0,
                'EX_taur_e': 0.0,
                'EX_tcb_e': 0.0,
                'EX_tcynt_e': 0.0, # Thiocynate (CNS)
                'EX_thm_e': 0.0,
                'EX_thr__L_e': 0.0,
                'EX_thym_e': 0.0,
                'EX_tmanur_e': 0.0,
                'EX_tnt_e': 0.0,
                'EX_tol_e': 0.0,
                'EX_tre_e': 0.0,
                'EX_trp__L_e': 0.0,
                'EX_tsul_e': 1000.0,
                'EX_ttdca_e': 0.0,
                'EX_ttdcea_e': 0.0,
                'EX_tyr__L_e': 0.0,
                'EX_ura_e': 0.0,
                'EX_urea_e': 0.0, # Urea
                'EX_uri_e': 0.0,
                'EX_vacc_e': 0.0,
                'EX_val__D_e': 0.0,
                'EX_val__L_e': 0.0,
                'EX_vanlt_e': 0.0,
                'EX_xan_e': 0.0,
                'EX_xtsn_e': 0.0,
                'EX_xyl__D_e': 0.0,
                'EX_zn2_e': 1000.0}

        model.medium = medium

        # Run the model
        sol = model.optimize()

        # Check that no biomass is produced
        self.assertEqual(sol.objective_value, 0)

if __name__ == '__main__':
    unittest.main()