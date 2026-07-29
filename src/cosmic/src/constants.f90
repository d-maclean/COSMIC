module constants

   integer, parameter:: dp = selected_real_kind(p = 15)

   real(dp), public, protected::&
      zsun = 0.014d0,&
      yr_to_day = 365.24d0,&
      au_to_rsun = 214.95d0,&
      yr_to_sec = 3.1557d7,&
      m_chandr = 1.44d0,&
      tiny = 1d-14,&
      acc1 = 3.920659d8,&
      k3 = 0.21d0,&
      mr23yr = 0.4311d0,&
      twopi = 2.d0*ACOS(-1.d0),&
      kw3 = 619.2d0,&
      wsun = 9.46d7,&
      wx = 9.46d8

   ! key evolutionary stages
   integer, public, protected::&
      initial_state = 1,&
      kstar_change = 2,&
      begin_RLOF = 3,&
      end_RLOF = 4,&
      contact = 5,&
      coalescence = 6,&
      begin_CE = 7,&
      end_CE = 8,&
      no_remnant_leftover = 8,&
      final_state = 10,&
      binary_disruption = 11,&
      begin_symbiotic_phase = 12,&
      end_symbiotic_phase = 13,&
      blue_straggler = 14,&
      supernova_primary = 15,&
      supernova_secondary = 16,&
      error_rlof_timeout = 100

   ! SSE star types
   ! code from Poojan Agrawal, METISSE
   integer, public, protected::&
      low_mass_MS = 0,&
      MS = 1,&
      HG = 2,&
      RGB = 3,&
      HeBurn = 4,&
      EAGB = 5,&
      TPAGB = 6,&
      He_MS = 7,&
      He_HG = 8,&
      He_GB = 9,&
      HeWD = 10,&
      CO_WD = 11,&
      ONeWD = 12,&
      NS = 13,&
      BH = 14,&
      Massless_REM = 15

   ! wind flag
   integer, public, protected::&
      ml_no_wind = -1,&
      ml_hurley_00 = 0,&
      ml_startrack_08 = 1,&
      ml_vink_01_05 = 2,&
      ml_vink_lbv = 3,&
      ml_vink_div3 = 5,&
      ml_bjorklund = 6,&
      ml_krticka_24 = 7

   ! lbvflag
   integer, public, protected::&
      lbv_off = 0,&
      lbv_hurley_00 = 1,&
      lbv_belczynski_08 = 2

   ! qcflag
   integer, public, protected::&
      qc_hurley_02 = 0,&
      qc_hurley_hjellming = 1,&
      qc_claeys_14 = 2,&
      qc_claeys_hjellming = 3,&
      qc_belczynski_08 = 4,&
      qc_neijsssel_20 = 5

   ! kickflag
   integer, public, protected::&
      kick_std = 1,&
      kick_giacobbo_20 = 2,&
      kick_giaccobo_20_ejecta_only = 3,&
      kick_bray_16 = 4,&
      kick_disberg_25 = 5,&
      kick_mandel_20 = 6
   ! negative values instead use Kiel & Hurley

   ! cekickflag
   integer, public, protected::&
      cekick_m_pre_sep_post = 0,&
      cekick_m_pre_sep_pre = 1,&
      cekick_m_post_sep_post = 2

   ! cehestarflag
   integer, public, protected::&
      ce_he_star_none = 0,&
      ce_he_star_period = 1,&
      ce_he_star_period_and_mass = 2

   !remnantflag
   integer, public, protected::&
      rem_hurley_00 = 0,&
      rem_belczynski_02 = 1,&
      rem_belczynski_08 = 2,&
      rem_fryer_12_rapid = 3,&
      rem_fryer_12_delayed = 4,&
      rem_mandel_20 = 5,&
      rem_maltsev_25 = 6,&
      rem_fryer_22 = 7

end module constants
