set(MODULES
  danilo_plot_combined
  danilo_datacards
  danilo_distributions
  danilo_check_scan
  danilo_SFapplication
  danilo_check_overlap
  danilo_plot2d_scan
  danilo_signal_scan
  danilo_cutflow_GGM
  danilo_acceptanceHist
  danilo_totsignalyield
  danilo_nDecays
  danilo_signal_scan_newTrig
  danilo_interpolateAcc
  danilo_plot_postfitBKG
  danilo_rValuediff
  danilo_GGM_scan
  danilo_plot1d_limit
  danilo_maxSensitivity
  danilo_plot_GGM
  danilo_plotGGMparameter
  danilo_plotOverlap
  danilo_plotSig
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
